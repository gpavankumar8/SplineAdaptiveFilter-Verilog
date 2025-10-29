`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 16-03-2023
// Module name: hsaflms_top.v
//////////////////////////////////////////////////////////////////////////////////

module hsaflms_top
#(parameter L_ORD = 32, Q_ORD = 4, WIDTH = 16, QP = 12, DelX_inv = 2, Q = 13)
(
    input clk,
    input reset, 
    input [WIDTH-1:0] signal_in,
    input [WIDTH-1:0] desired_in,
    output [WIDTH-1:0] filter_out_d,
    output reg [WIDTH-1:0] error_d
);

    localparam WRET = 0, NON_PIP = 2, E_PIP = 0;
    localparam LIN_PIP_FIR = 1 + (($clog2(L_ORD)-1)/4);
    localparam NON_PIP_FIR = 1 + (($clog2(Q_ORD)-1)/4);
    localparam DOTP_PIP = 1 + (($clog2(L_ORD)-1)/4);
    localparam RET = NON_PIP + NON_PIP_FIR + LIN_PIP_FIR + E_PIP + WRET;
    localparam W_PIP = LIN_PIP_FIR + E_PIP;
    localparam Q_PIP = NON_PIP_FIR + LIN_PIP_FIR + E_PIP - DOTP_PIP;
    localparam LIN_ADD_LEN = 2**$clog2(L_ORD);
    localparam LIN_ADD_ZERO_PAD = LIN_ADD_LEN - L_ORD ;
    localparam NONLIN_ADD_LEN = 2**$clog2(Q_ORD);
    localparam NONLIN_ADD_ZERO_PAD = NONLIN_ADD_LEN - Q_ORD;

    reg signed  [WIDTH-1:0] signal_in_d, desired_in_d;
    wire signed [WIDTH-1:0] lin_filter_in[L_ORD+W_PIP-1:0], weight[L_ORD-1:0], mu_w_error, mu_w_error_d, error_w_rnd, mu_q_error, mu_q_error_d, error, error_q_rnd, adder_out, filter_out;
    wire signed [WIDTH-1:0] adder_tree_out, adder_tree_out_rnd, desired_in_ret_d;

    wire signed [L_ORD*WIDTH-1:0]     tap_out_packed, tap_out_packed_d, weight_packed, weight_packed_vmm;
    wire signed [(L_ORD+W_PIP-1)*WIDTH-1:0] lin_filter_in_packed;

    wire signed [(L_ORD-1)*WIDTH-1:0] nonl_filter_in_packed[Q_ORD-1:0];
    wire signed [Q_ORD*WIDTH-1:0] a_out_packed, a_weight_packed;
    wire signed [WIDTH-1:0] nonl_filter_in[Q_ORD-1:0][L_ORD-1:0], s_out, a_weight[Q_ORD-1:0];
    wire signed [WIDTH-1:0] weight_pip_prod[Q_ORD-1:0];

    // Nonlinear spline interpolation module

    localparam SPAN_WIDTH = $clog2(Q+Q_ORD);
    
    // wire signed [WIDTH-1:0] u, u_d, u_2d, u2, u2_d, u2_2d, u3, u3_d, x_div_DelX, u_vec_C[Q_ORD-1:0], u_vec_C_d[Q_ORD-1:0], q_weight[Q_ORD-1:0], q_weight_old[Q_ORD-1:0], q_update[Q_ORD-1:0];
    // wire signed [Q_ORD*WIDTH-1:0] u_vec_C_packed, u_vec_C_out_packed, q_weight_packed, q_weight_old_packed, q_update_packed, s_out_packed, s_out_packed_d;
    // wire signed [(L_ORD+Q_PIP-2)*WIDTH-1:0] u_vec_C_pipeline[Q_ORD-1:0];
    // wire signed [L_ORD*WIDTH-1:0] u_vec_C_pipeline_span_cmp[Q_ORD-1:0], u_vec_C_pipeline_span_cmp_d[Q_ORD-1:0], span_ind_cmp_expand;
    // wire signed [2*WIDTH-1:0] u2_full, u3_full, u2_rnd, u3_rnd;
    // wire [SPAN_WIDTH-1:0] span_ind, span_ind_read, span_ind_write, span_ind_write_d;
    // wire [(L_ORD-1)*SPAN_WIDTH-1:0] span_ind_pipeline;
    // wire span_ind_cmp[L_ORD-1:0];

    // // Generate u and span index (j)
    // assign x_div_DelX = signal_in_d <<< DelX_inv;
    // assign u = x_div_DelX[QP-1:0];
    // assign span_ind = $signed(x_div_DelX[WIDTH-1:QP]) + $signed((Q-1)/2);
    // DelayNUnit #(WIDTH, 1) U_PIP(clk, reset, u, u_d);
    // DelayNUnit #(WIDTH, 1) U_D_PIP(clk, reset, u_d, u_2d);

    // // Generate u_vec_C vector
    // assign u2_full = u  * u;
    // assign u2_rnd = u2_full + (1'b1 << (QP-1));
    // assign u2 = u2_rnd[QP+:WIDTH];
    // DelayNUnit #(WIDTH, 1) U2_PIP(clk, reset, u2, u2_d);
    // DelayNUnit #(WIDTH, 1) U2_D_PIP(clk, reset, u2_d, u2_2d);

    // assign u3_full = u2_d * u_d;
    // assign u3_rnd = u3_full + (1'b1 << (QP-1));
    // assign u3 = u3_rnd[QP+:WIDTH];
    // DelayNUnit #(WIDTH, 1) U3_PIP(clk, reset, u3, u3_d);

    // // Hardcoded u_vec_C for CR spline matrix and P = 3 (Q_ORD = 4)
    // // assign u_vec_C[0] = (-1*u3 + 2*u2 - u) >>> 1;
    // // assign u_vec_C[1] = ( 3*u3 - 5*u2 + (2'b10 << QP)) >>> 1;
    // // assign u_vec_C[2] = (-3*u3 + 4*u2 + u) >>> 1;
    // // assign u_vec_C[3] = (   u3 -   u2    ) >>> 1;

    // assign u_vec_C[0] = (-1*u3_d + 2*u2_2d - u_2d) >>> 1;
    // assign u_vec_C[1] = ( 3*u3_d - 5*u2_2d + (2'b10 << QP)) >>> 1;
    // assign u_vec_C[2] = (-3*u3_d + 4*u2_2d + u_2d) >>> 1;
    // assign u_vec_C[3] = (   u3_d -   u2_2d       ) >>> 1;

    wire signed [WIDTH-1:0] u_vec_C[Q_ORD-1:0], u_vec_C_d[Q_ORD-1:0], q_weight[Q_ORD-1:0], q_weight_old[Q_ORD-1:0], q_update[Q_ORD-1:0];
    wire signed [Q_ORD*WIDTH-1:0] u_vec_C_packed, u_vec_C_out_packed, q_weight_packed, q_weight_old_packed, q_update_packed, s_out_packed, s_out_packed_d;
    wire signed [(L_ORD+Q_PIP-2)*WIDTH-1:0] u_vec_C_pipeline[Q_ORD-1:0];
    wire signed [L_ORD*WIDTH-1:0] u_vec_C_pipeline_span_cmp[Q_ORD-1:0], u_vec_C_pipeline_span_cmp_d[Q_ORD-1:0], span_ind_cmp_expand, span_ind_cmp_expand_d;
    wire [SPAN_WIDTH-1:0] span_ind, span_ind_read, span_ind_write, span_ind_write_d;
    wire [(L_ORD-1)*SPAN_WIDTH-1:0] span_ind_pipeline;
    wire span_ind_cmp[L_ORD-1:0], span_ind_cmp_2d[L_ORD-1:0];
    wire [L_ORD-1:0] span_ind_cmp_d;

    u_vec_span_gen #(WIDTH, Q_ORD, QP, Q, DelX_inv, SPAN_WIDTH) U_VEC_GEN(clk, reset, signal_in_d, u_vec_C_out_packed, span_ind);
    
    // Store u_vec_C in pipeline
    genvar f;
    generate
        for (f = 0; f < Q_ORD; f = f+1)
        begin:pipe
            DelayNUnit #(WIDTH, 1) U_VEC_C_DEL(clk, reset, u_vec_C[f], u_vec_C_d[f]);               // Retiming delay after U_vec * C operation
            pipeline #(WIDTH, (L_ORD+Q_PIP-2)) U_VEC_C_PIP(clk, reset, u_vec_C_d[f], u_vec_C_pipeline[f]);
        end
    endgenerate

    // Q weight selection and update
    // assign span_ind_read = span_ind;
    // assign span_ind_write = span_ind;
    
    DelayNUnit #(SPAN_WIDTH, NON_PIP-1) SPN_READ_DEL(clk, reset, span_ind, span_ind_read);               // Retiming delay for span_ind_read to q_weight_controller
    // DelayNUnit #(SPAN_WIDTH, Q_PIP-1) SPN_WRITE_DEL(clk, reset, span_ind_read , span_ind_write);               // Retiming delay for span_ind_write to q_weight_controller
    assign span_ind_write = span_ind_read;
    DelayNUnit #(SPAN_WIDTH, DOTP_PIP+2) SPN_WRITE_DEL2(clk, reset, span_ind_write, span_ind_write_d);

    pipeline #(SPAN_WIDTH, (L_ORD-1)) SPAN_IND_PIP(clk, reset, span_ind_write, span_ind_pipeline);
    q_weight_controller_d #(WIDTH, Q, Q_ORD) Q_WEIGHTS(.clk(clk), .reset(reset), .span_ind_write_d(span_ind_write_d), .q_update_packed(q_update_packed), .span_ind_write(span_ind_write_d), .q_weight_old_packed(q_weight_old_packed), .span_ind_read(span_ind_read), .q_weight_packed_out(q_weight_packed));

    // Nonlinear filter output
    fir_taps #(WIDTH, QP, Q_ORD) SPLINE_FIR(u_vec_C_packed, q_weight_packed, s_out_packed); 
    // DelayNUnit #(Q_ORD*WIDTH, 1) NONL_MULT_OUT_PIP(clk, reset, s_out_packed, s_out_packed_d);              // Retiming delay after NONL FIR multipliers
    // adder_tree #(WIDTH, Q_ORD) SPLINE_ADD(.adder_tree_in_packed(s_out_packed_d), .adder_tree_out(s_out));
    adder_tree #(WIDTH, NONLIN_ADD_LEN) SPLINE_ADD(.clk(clk), .reset(reset), .adder_tree_in_packed({s_out_packed,{(NONLIN_ADD_ZERO_PAD*WIDTH){1'b0}}}), .adder_tree_out(s_out)); 
    DelayNUnit #(WIDTH, 1) NONL_FILT_OUT_PIP(clk, reset, s_out, lin_filter_in[0]);              // Retiming delay after NONL adder tree

    //assign lin_filter_in[0] = s_out;

    // Linear filter
    pipeline #(WIDTH, (L_ORD+W_PIP-1)) S_PIP(clk, reset, lin_filter_in[0], lin_filter_in_packed);
    fir_taps #(WIDTH, QP, L_ORD) LIN_FIR({lin_filter_in_packed[(L_ORD-1)*WIDTH-1:0], lin_filter_in[0]}, weight_packed, tap_out_packed); 
    // fir_taps #(WIDTH, QP, 1) FIR_CONST((1'b1 << QP), weight_const, tap_out_const);
    // DelayNUnit #(L_ORD*WIDTH, 1) MULT_PIP(clk, reset, tap_out_packed, tap_out_packed_d);      // Retiming delay after multipliers

    // adder_tree #(WIDTH, L_ORD) LIN_ADD(.adder_tree_in_packed(tap_out_packed), .adder_tree_out(adder_out)); 
    adder_tree #(WIDTH, LIN_ADD_LEN) LIN_ADD(.clk(clk), .reset(reset), .adder_tree_in_packed({tap_out_packed,{(LIN_ADD_ZERO_PAD*WIDTH){1'b0}}}), .adder_tree_out(adder_out)); 
    assign filter_out = adder_out;
    DelayNUnit #(WIDTH, 1) FILT_OUT_PIP(clk, reset, filter_out, filter_out_d);

    // Compute error
    error_compute #(WIDTH) EC( .desired_in(desired_in_ret_d), .filter_out(filter_out_d), .error(error));

    // Linear Weight update
    assign error_w_rnd = error + (1<<(7-1));
    assign mu_w_error = error_w_rnd >>> 7;       // Multiply by mu = 0.0078125 (1/(2^7));

    // DelayNUnit #(WIDTH, 1) ERR_PIP_W(clk, reset, mu_w_error, mu_w_error_d);              // Retiming delay after error_w

    assign error_q_rnd = error + (1<<(7-1));
    assign mu_q_error = error_q_rnd >>> 7;       // Multiply by mu = 0.015625 (1/(2^6));

    // DelayNUnit #(WIDTH, 1) ERR_PIP_Q(clk, reset, mu_q_error, mu_q_error_d);              // Retiming delay after error_q

    genvar w;
    generate            
        for ( w = 0; w < L_ORD; w = w+1 ) 
        begin: weights_l
            w_update_d #(WIDTH, QP) WUB( clk, reset, mu_w_error, lin_filter_in[w+W_PIP], weight[w]); 
        end
    endgenerate

    // Weight for constant input
    //w_update #(WIDTH, QP) WUB_CONST( clk, reset, mu_w_error, (1'b1 << QP), weight_const); 

    // Nonlinear (a) weight update

    // dot_product #(WIDTH, QP, L_ORD) WEIGHT0_PIP({nonl_filter_in_packed[0], nonl_filter_in[0][0]}, weight_packed, weight_pip_prod[0]);
    // w_update #(WIDTH, QP, 16'h0001) A_WUB0( clk, reset, mu_a_error, weight_pip_prod[0], a_weight[0]); 

    assign span_ind_cmp_d[0]  = 1'b1;
    assign span_ind_cmp_expand[WIDTH-1:0] = {WIDTH{span_ind_cmp_d[0]}};

    genvar a;
    generate            
        for ( a = 1; a < L_ORD; a = a+1 ) 
        begin: comparators
            assign span_ind_cmp[a] = (span_ind_write == span_ind_pipeline[((a-1)*SPAN_WIDTH)+:SPAN_WIDTH]);
            DelayNUnit #(1, 1) SPAN_CMP_PIP(clk, reset, span_ind_cmp[a], span_ind_cmp_d[a]);
            assign span_ind_cmp_expand[WIDTH*a+:WIDTH] = {WIDTH{span_ind_cmp_d[a]}};
        end
    endgenerate

    DelayNUnit #(L_ORD*WIDTH, 1) SPAN_CMP_EXP_PIP(clk, reset, span_ind_cmp_expand, span_ind_cmp_expand_d);

    // Weight input to VMM stage delays
    // DelayNUnit #(L_ORD*WIDTH, 1+LIN_PIP_FIR-DOTP_PIP) WEIGHT_VMM_PIP(clk, reset, weight_packed, weight_packed_vmm);
    DelayUnit_CG #(L_ORD, WIDTH) WEIGHT_VMM_PIP(clk, reset, span_ind_cmp_d, weight_packed, weight_packed_vmm);

    generate            
        for ( a = 0; a < Q_ORD; a = a+1 ) 
        begin: weights_nonl
            // assign u_vec_C_pipeline_span_cmp[a] = {u_vec_C_pipeline[a], u_vec_C_d[a]} & span_ind_cmp_expand;
            assign u_vec_C_pipeline_span_cmp[a] = {u_vec_C_pipeline[a], u_vec_C_d[a]};
            DelayUnit_CG #(L_ORD, WIDTH) ERR_PIP_W(clk, reset, span_ind_cmp_d, u_vec_C_pipeline_span_cmp[a], u_vec_C_pipeline_span_cmp_d[a]);
            // DelayNUnit #(L_ORD*WIDTH, 1) ERR_PIP_W(clk, reset, u_vec_C_pipeline_span_cmp[a], u_vec_C_pipeline_span_cmp_d[a]);
            // dot_product #(WIDTH, QP, L_ORD) WEIGHT_PIP({u_vec_C_pipeline_span_cmp[a], u_vec_C[a]}, weight_packed, weight_pip_prod[a]);
            dot_product_d_gated #(WIDTH, QP, L_ORD) WEIGHT_PIP(clk, reset, span_ind_cmp_expand_d, u_vec_C_pipeline_span_cmp_d[a], weight_packed_vmm, weight_pip_prod[a]);
            w_update_term_d #(WIDTH, QP) A_WUB( clk, reset, mu_q_error, weight_pip_prod[a], q_weight_old[a], q_update[a]); 
        end
    endgenerate

    always @ (posedge clk)
    begin
        if (reset)
        begin
            signal_in_d <= 0;
            desired_in_d <= 0;           
        end
        else
        begin
            signal_in_d <= signal_in;
            desired_in_d <= desired_in;
        end
    end

    always @ (posedge clk)
    begin
        if (reset)
        begin
            error_d <= 0;
        end
        else
        begin
            error_d <= error;
        end
    end
 
    // Delay desired signal based on retiming delays in signal_in path
    DelayNUnit #(WIDTH, RET-WRET-E_PIP) DES_PIP(clk, reset, desired_in_d, desired_in_ret_d);

    // 2D and 1D array conversions
    // Unpack nonl_filter_in pipeline
    genvar ind_l, ind_q;
    generate
        for ( ind_q = 0; ind_q < Q_ORD; ind_q=ind_q+1 )
        begin:filter_pack
            for ( ind_l = 0; ind_l < L_ORD-1; ind_l=ind_l+1 ) 
            begin:nonl_filter_inner
                assign nonl_filter_in[ind_q][ind_l+1] = nonl_filter_in_packed[ind_q][WIDTH*ind_l+:WIDTH];
            end
        end
    endgenerate

    genvar ind;
    generate
        for ( ind = 0; ind < L_ORD; ind=ind+1 )
        begin:weight_pack
            assign weight_packed[WIDTH*ind+:WIDTH] = weight[ind];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:q_weight_unpack
            assign q_weight[ind] = q_weight_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:q_weight_old_unpack
            assign q_weight_old[ind] = q_weight_old_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:q_update_pack
            assign q_update_packed[WIDTH*ind+:WIDTH] = q_update[ind];
        end
    endgenerate

    generate
        for ( ind = 0; ind < L_ORD+W_PIP-1; ind=ind+1 )
        begin:lin_pack
            assign lin_filter_in[ind+1] = lin_filter_in_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:u_vec_C_pack
            assign u_vec_C_packed[WIDTH*ind+:WIDTH] = u_vec_C_d[ind];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:u_vec_C_unpack
            assign u_vec_C[ind] = u_vec_C_out_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate


endmodule