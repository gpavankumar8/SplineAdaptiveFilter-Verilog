`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 03-03-2023
// Module name: hsaflms_top.v
//////////////////////////////////////////////////////////////////////////////////

module hsaflms_top
#(parameter L_ORD = 32, Q_ORD = 4, WIDTH = 16, QP = 12, DelX_inv = 2, Q = 13)
(
    input clk,
    input reset, 
    input [WIDTH-1:0] signal_in,
    input [WIDTH-1:0] desired_in,
    output [WIDTH-1:0] filter_out,
    output [WIDTH-1:0] error
);

    reg signed  [WIDTH-1:0] signal_in_d, desired_in_d;
    wire signed [WIDTH-1:0] lin_filter_in[L_ORD-1:0], weight[L_ORD-1:0], mu_w_error, error_w_rnd, mu_q_error, error_q_rnd, adder_out;
    wire signed [WIDTH-1:0] adder_tree_out, adder_tree_out_rnd;

    wire signed [L_ORD*WIDTH-1:0]     tap_out_packed, weight_packed;
    wire signed [(L_ORD-1)*WIDTH-1:0] lin_filter_in_packed;

    wire signed [Q_ORD*WIDTH-1:0] nonl_x_out_packed;
    wire signed [WIDTH-1:0]       nonl_x_out[Q_ORD-1:0];

    wire signed [WIDTH-1:0] weight_const, tap_out_const;

    wire signed [(L_ORD-1)*WIDTH-1:0] nonl_filter_in_packed[Q_ORD-1:0];
    wire signed [Q_ORD*WIDTH-1:0] a_out_packed, a_weight_packed;
    wire signed [WIDTH-1:0] nonl_filter_in[Q_ORD-1:0][L_ORD-1:0], s_out, a_weight[Q_ORD-1:0];
    wire signed [WIDTH-1:0] weight_pip_prod[Q_ORD-1:0];

    // // Nonlinear Phi mapping
    // nonl_PhiMap #(Q_ORD, WIDTH, QP) Phi(.x_in(signal_in_d), .nonl_x_out_packed(nonl_x_out_packed));

    // // Filter pipeline
    // genvar f;
    // generate
    //     for (f = 0; f < Q_ORD; f = f+1)
    //     begin:pipe
    //         assign nonl_filter_in[f][0] = nonl_x_out[f];
    //         pipeline #(WIDTH, (L_ORD-1)) NONL_PIP(clk, reset, nonl_filter_in[f][0], nonl_filter_in_packed[f]);
    //     end
    // endgenerate

    // // Non-linearity weights (a)
    // fir_taps #(WIDTH, QP, Q_ORD) NONL_FIR(nonl_x_out_packed, a_weight_packed, a_out_packed); 
    // adder_tree #(WIDTH, Q_ORD+1) NONL_ADD1(.adder_tree_in_packed({a_out_packed,{WIDTH{1'b0}}}), .adder_tree_out(s_out)); 

    // Nonlinear spline interpolation module

    localparam SPAN_WIDTH = $clog2(Q+Q_ORD);
    
    wire signed [WIDTH-1:0] u, u2, u3, x_div_DelX, u_vec_C[Q_ORD-1:0], q_weight[Q_ORD-1:0], q_update[Q_ORD-1:0];
    wire signed [Q_ORD*WIDTH-1:0] u_vec_C_packed, q_weight_packed, q_update_packed, s_out_packed;
    wire signed [(L_ORD-1)*WIDTH-1:0] u_vec_C_pipeline[Q_ORD-1:0], u_vec_C_pipeline_span_cmp[Q_ORD-1:0], span_ind_cmp_expand;
    wire signed [2*WIDTH-1:0] u2_full, u3_full, u2_rnd, u3_rnd;
    wire [SPAN_WIDTH-1:0] span_ind, span_ind_read, span_ind_write;
    wire [(L_ORD-1)*SPAN_WIDTH-1:0] span_ind_pipeline;
    wire span_ind_cmp[(L_ORD-1)-1:0];

    // Generate u and span index (j)
    assign x_div_DelX = signal_in_d <<< DelX_inv;
    assign u = x_div_DelX[QP-1:0];
    assign span_ind = $signed(x_div_DelX[WIDTH-1:QP]) + $signed((Q-1)/2);

    // Generate u_vec_C vector
    assign u2_full = u  * u;
    assign u2_rnd = u2_full + (1'b1 << (QP-1));
    assign u2 = u2_rnd[QP+:WIDTH];

    assign u3_full = u2 * u;
    assign u3_rnd = u3_full + (1'b1 << (QP-1));
    assign u3 = u3_rnd[QP+:WIDTH];

    // Hardcoded u_vec_C for CR spline matrix and P = 3 (Q_ORD = 4)
    assign u_vec_C[0] = (-1*u3 + 2*u2 - u) >>> 1;
    assign u_vec_C[1] = ( 3*u3 - 5*u2 + (2'b10 << QP)) >>> 1;
    assign u_vec_C[2] = (-3*u3 + 4*u2 + u) >>> 1;
    assign u_vec_C[3] = (   u3 -   u2    ) >>> 1;

    // Store u_vec_C in pipeline
    genvar f;
    generate
        for (f = 0; f < Q_ORD; f = f+1)
        begin:pipe
            pipeline #(WIDTH, (L_ORD-1)) U_VEC_C_PIP(clk, reset, u_vec_C[f], u_vec_C_pipeline[f]);
        end
    endgenerate

    // Q weight selection and update
    assign span_ind_read = span_ind;
    assign span_ind_write = span_ind;
    pipeline #(SPAN_WIDTH, (L_ORD-1)) SPAN_IND_PIP(clk, reset, span_ind, span_ind_pipeline);
    q_weight_controller #(WIDTH, Q, Q_ORD) Q_WEIGHTS(.clk(clk), .reset(reset), .span_ind_write(span_ind_write), .q_update_packed(q_update_packed), .span_ind_read(span_ind_read), .q_weight_packed_out(q_weight_packed));

    // Nonlinear filter output
    fir_taps #(WIDTH, QP, Q_ORD) SPLINE_FIR(u_vec_C_packed, q_weight_packed, s_out_packed); 
    adder_tree #(WIDTH, Q_ORD) SPLINE_ADD(.adder_tree_in_packed(s_out_packed), .adder_tree_out(s_out));

    assign lin_filter_in[0] = s_out;

    // Linear filter
    pipeline #(WIDTH, (L_ORD-1)) S_PIP(clk, reset, lin_filter_in[0], lin_filter_in_packed);
    fir_taps #(WIDTH, QP, L_ORD) LIN_FIR({lin_filter_in_packed, lin_filter_in[0]}, weight_packed, tap_out_packed); 
    // fir_taps #(WIDTH, QP, 1) FIR_CONST((1'b1 << QP), weight_const, tap_out_const);

    adder_tree #(WIDTH, L_ORD) LIN_ADD(.adder_tree_in_packed(tap_out_packed), .adder_tree_out(adder_out)); 
    assign filter_out = adder_out;

    // Compute error
    error_compute #(WIDTH) EC( .desired_in(desired_in_d), .filter_out(filter_out), .error(error));

    // Linear Weight update
    assign error_w_rnd = error + (1<<(7-1));
    assign mu_w_error = error_w_rnd >>> 7;       // Multiply by mu = 0.0078125 (1/(2^7));

    assign error_q_rnd = error + (1<<(7-1));
    assign mu_q_error = error_q_rnd >>> 7;       // Multiply by mu = 0.015625 (1/(2^6));

    genvar w;
    generate            
        for ( w = 0; w < L_ORD; w = w+1 ) 
        begin: weights_l
            w_update #(WIDTH, QP) WUB( clk, reset, mu_w_error, lin_filter_in[w], weight[w]); 
        end
    endgenerate

    // Weight for constant input
    //w_update #(WIDTH, QP) WUB_CONST( clk, reset, mu_w_error, (1'b1 << QP), weight_const); 

    // Nonlinear (a) weight update

    // dot_product #(WIDTH, QP, L_ORD) WEIGHT0_PIP({nonl_filter_in_packed[0], nonl_filter_in[0][0]}, weight_packed, weight_pip_prod[0]);
    // w_update #(WIDTH, QP, 16'h0001) A_WUB0( clk, reset, mu_a_error, weight_pip_prod[0], a_weight[0]); 

    genvar a;
    generate            
        for ( a = 0; a < L_ORD-1; a = a+1 ) 
        begin: comparators
            assign span_ind_cmp[a] = (span_ind_write == span_ind_pipeline[SPAN_WIDTH*a+:SPAN_WIDTH]);
            assign span_ind_cmp_expand[WIDTH*a+:WIDTH] = {WIDTH{span_ind_cmp[a]}};
        end
    endgenerate

    generate            
        for ( a = 0; a < Q_ORD; a = a+1 ) 
        begin: weights_nonl
            assign u_vec_C_pipeline_span_cmp[a] = u_vec_C_pipeline[a] & span_ind_cmp_expand;
            dot_product #(WIDTH, QP, L_ORD) WEIGHT_PIP({u_vec_C_pipeline_span_cmp[a], u_vec_C[a]}, weight_packed, weight_pip_prod[a]);
            w_update_term #(WIDTH, QP) A_WUB( clk, reset, mu_q_error, weight_pip_prod[a], q_weight[a], q_update[a]); 
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
        begin:q_update_pack
            assign q_update_packed[WIDTH*ind+:WIDTH] = q_update[ind];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:nonl_out_pack
            assign nonl_x_out[ind] = nonl_x_out_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < L_ORD-1; ind=ind+1 )
        begin:lin_pack
            assign lin_filter_in[ind+1] = lin_filter_in_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:u_vec_C_pack
            assign u_vec_C_packed[WIDTH*ind+:WIDTH] = u_vec_C[ind];
        end
    endgenerate

endmodule