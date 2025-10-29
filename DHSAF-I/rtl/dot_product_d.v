`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 17-10-2022
// Module name: dot_product_d.v
//////////////////////////////////////////////////////////////////////////////////

module dot_product_d
#(parameter WIDTH = 16, QP = 12, LEN = 8)
(
    input clk,
    input reset,
    input [LEN*WIDTH-1:0] vec1_packed,
    input [LEN*WIDTH-1:0] vec2_packed,
    output [WIDTH-1:0] dotp_out
);

    wire signed [LEN*WIDTH-1:0] mult_out_packed, mult_out_packed_d;
    wire signed [WIDTH-1:0]     dotp;

    fir_taps   #(WIDTH, QP, LEN) DOTP(vec1_packed, vec2_packed, mult_out_packed);
    // DelayNUnit #(LEN*WIDTH, 1) MULT_PIP(clk, reset, mult_out_packed, mult_out_packed_d);              // Retiming delay after multipliers

    adder_tree #(WIDTH, LEN) DOT_ADD(.clk(clk), .reset(reset), .adder_tree_in_packed(mult_out_packed), .adder_tree_out(dotp));
    DelayNUnit #(WIDTH, 1) ADD_PIP(clk, reset, dotp, dotp_out);                                       // Retiming delay after adder tree

endmodule