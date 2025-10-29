`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 11-08-2022
// Module name: error_compute.v
//////////////////////////////////////////////////////////////////////////////////

module error_compute
#(parameter WIDTH = 16)
(
    input [WIDTH-1:0] desired_in,
    input [WIDTH-1:0] filter_out,
    output [WIDTH-1:0] error
);

    assign error = desired_in - filter_out;

endmodule