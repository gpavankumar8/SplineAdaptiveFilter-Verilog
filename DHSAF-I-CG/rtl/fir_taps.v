`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 22-08-2022
// Module name: fir_taps.v
//////////////////////////////////////////////////////////////////////////////////

module fir_taps
#(parameter WIDTH = 16, QP = 12, ORD = 64)
(
    input [ORD*WIDTH-1:0] filter_in_packed,
    input [ORD*WIDTH-1:0] weight_in_packed,
    output [ORD*WIDTH-1:0] tap_out_packed
);

    wire signed [2*WIDTH-1:0] tap_out_full[ORD-1:0], tap_out_rnd[ORD-1:0];
    wire signed [WIDTH-1:0] filter_in[ORD-1:0], weight[ORD-1:0], tap_out[ORD-1:0];

    genvar i;

    generate
        for(i = 0; i < ORD; i = i + 1)
        begin:taps
            assign tap_out_full[i] = $signed(filter_in[i]) * $signed(weight[i]);
            assign tap_out_rnd[i] = tap_out_full[i] + (1'b1 << (QP-1));
            assign tap_out[i] = tap_out_rnd[i][QP+:WIDTH];
            // tap_multiply #(WIDTH, QP) TAP(filter_in[i], weight[i], tap_out[i]);
        end
    endgenerate


    genvar ind;
    generate
        for ( ind = 0; ind < ORD; ind=ind+1 )
        begin:filter_pack
            assign filter_in[ind] = filter_in_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < ORD; ind=ind+1 )
        begin:weight_pack
            assign weight[ind] = weight_in_packed[WIDTH*ind+:WIDTH];
        end
    endgenerate

    generate
        for ( ind = 0; ind < ORD; ind=ind+1 )
        begin:tap_pack
            assign tap_out_packed[WIDTH*ind+:WIDTH] = tap_out[ind];
        end
    endgenerate

endmodule