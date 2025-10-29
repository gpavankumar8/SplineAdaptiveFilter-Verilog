`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 18-03-2023
// Module name: q_weight_controller_d.v
//////////////////////////////////////////////////////////////////////////////////

module q_weight_controller_d
#(parameter WIDTH = 16, Q = 13, Q_ORD = 4)
(
    input clk,
    input reset,
    input [$clog2(Q+Q_ORD)-1:0] span_ind_write_d,
    input [Q_ORD*WIDTH-1:0] q_update_packed,
    input [$clog2(Q+Q_ORD)-1:0] span_ind_write,
    output [Q_ORD*WIDTH-1:0] q_weight_old_packed,
    input [$clog2(Q+Q_ORD)-1:0] span_ind_read,
    output [Q_ORD*WIDTH-1:0] q_weight_packed_out
);

    reg [WIDTH-1:0] q_weight[Q+Q_ORD-1:0];
    wire [(Q+Q_ORD)*WIDTH-1:0] q_weight_packed;
    reg [(Q+Q_ORD)*WIDTH-1:0] q_weight_mux_out_packed;
    wire [WIDTH-1:0] q_weight_mux_out[Q+Q_ORD-1:0];
    wire [Q_ORD*WIDTH-1:0] q_weight_packed_out_mux;

    // Initialize and write Q weights

    always @ (posedge clk)
    begin
        if(reset)
        begin
            q_weight[0] <= 16'hE800;
            q_weight[1] <= 16'hEC00;
            q_weight[2] <= 16'hF000;
            q_weight[3] <= 16'hF400;
            q_weight[4] <= 16'hF800;
            q_weight[5] <= 16'hFC00;
            q_weight[6] <= 16'h0000;
            q_weight[7] <= 16'h0400;
            q_weight[8] <= 16'h0800;
            q_weight[9] <= 16'h0C00;
            q_weight[10] <= 16'h1000;
            q_weight[11] <= 16'h1400;
            q_weight[12] <= 16'h1800;
            q_weight[13] <= 16'h1C00;
            q_weight[14] <= 16'h2000;
            q_weight[15] <= 16'h2400;
            q_weight[16] <= 16'h2800;
        end
        else
        begin
            q_weight[0] <= q_weight_mux_out[0];
            q_weight[1] <= q_weight_mux_out[1];
            q_weight[2] <= q_weight_mux_out[2];
            q_weight[3] <= q_weight_mux_out[3];
            q_weight[4] <= q_weight_mux_out[4];
            q_weight[5] <= q_weight_mux_out[5];
            q_weight[6] <= q_weight_mux_out[6];
            q_weight[7] <= q_weight_mux_out[7];
            q_weight[8] <= q_weight_mux_out[8];
            q_weight[9] <= q_weight_mux_out[9];
            q_weight[10] <= q_weight_mux_out[10];
            q_weight[11] <= q_weight_mux_out[11];
            q_weight[12] <= q_weight_mux_out[12];
            q_weight[13] <= q_weight_mux_out[13];
            q_weight[14] <= q_weight_mux_out[14];
            q_weight[15] <= q_weight_mux_out[15];
            q_weight[16] <= q_weight_mux_out[16];
        end
    end

    // genvar i;
    // generate
    //     for( i = 0; i < (Q+Q_ORD); i = i+1 )
    //     begin:stage
    //         always @ (posedge clk)
    //             if( !reset )
    //                 q_weight[i] <= q_weight_mux_out[i];            
    //    end
    // endgenerate

    always @ (*)
    begin
        case(span_ind_write_d)
            5'd00: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD-0)*WIDTH],q_update_packed};
            5'd01: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+1)*WIDTH],q_update_packed, q_weight_packed[1*WIDTH-1:0]};
            5'd02: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+2)*WIDTH],q_update_packed, q_weight_packed[2*WIDTH-1:0]};
            5'd03: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+3)*WIDTH],q_update_packed, q_weight_packed[3*WIDTH-1:0]};
            5'd04: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+4)*WIDTH],q_update_packed, q_weight_packed[4*WIDTH-1:0]};
            5'd05: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+5)*WIDTH],q_update_packed, q_weight_packed[5*WIDTH-1:0]};
            5'd06: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+6)*WIDTH],q_update_packed, q_weight_packed[6*WIDTH-1:0]};
            5'd07: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+7)*WIDTH],q_update_packed, q_weight_packed[7*WIDTH-1:0]};
            5'd08: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+8)*WIDTH],q_update_packed, q_weight_packed[8*WIDTH-1:0]};
            5'd09: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+9)*WIDTH],q_update_packed, q_weight_packed[9*WIDTH-1:0]};
            5'd10: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+10)*WIDTH],q_update_packed, q_weight_packed[10*WIDTH-1:0]};
            5'd11: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+11)*WIDTH],q_update_packed, q_weight_packed[11*WIDTH-1:0]};
            5'd12: q_weight_mux_out_packed = {q_weight_packed[(Q+Q_ORD)*WIDTH-1:(Q_ORD+12)*WIDTH],q_update_packed, q_weight_packed[12*WIDTH-1:0]};
            5'd13: q_weight_mux_out_packed = {q_update_packed, q_weight_packed[13*WIDTH-1:0]};
            default: q_weight_mux_out_packed = {((Q+Q_ORD)*WIDTH-1){1'bx}};
        endcase
    end

    // Read Q weights
    assign q_weight_packed_out_mux = {q_weight[span_ind_read+3], q_weight[span_ind_read+2], q_weight[span_ind_read+1], q_weight[span_ind_read]};
    DelayNUnit #(Q_ORD*WIDTH, 1) Q_OUT_DEL(clk, reset, q_weight_packed_out_mux, q_weight_packed_out);               // Registered output

    // Read Q weights for weight update
    assign q_weight_old_packed = {q_weight[span_ind_write+3], q_weight[span_ind_write+2], q_weight[span_ind_write+1], q_weight[span_ind_write]};


    // Packing and unpacking arrays

    genvar ind;
    generate
        for ( ind = 0; ind < (Q+Q_ORD); ind=ind+1 )
        begin:q_weight_pack
            assign q_weight_packed[WIDTH*ind+:WIDTH] = q_weight[ind];
        end
    endgenerate

    generate
        for ( ind = 0; ind < (Q+Q_ORD); ind=ind+1 )
        begin:q_weight_mux_unpack
            assign q_weight_mux_out[ind] = q_weight_mux_out_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

endmodule