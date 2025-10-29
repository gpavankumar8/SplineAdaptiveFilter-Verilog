`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 11-08-2022
// Module name: w_update_d.v
//////////////////////////////////////////////////////////////////////////////////

module w_update_d
#(parameter WIDTH = 16, QP = 12, RESET_VAL = {WIDTH{1'b0}})
(
    input clk,
    input reset,
    input [WIDTH-1:0] mu_error,
    input [WIDTH-1:0] x_n,
    output reg [WIDTH-1:0] weight
);

    wire signed [2*WIDTH-1:0] x_n_error_full, x_n_error_rnd; 
    wire signed [WIDTH-1:0] x_n_error, new_weight, x_n_error_d;
    
    //DelayNUnit #(WIDTH, 1) IN_PIP(clk, reset, mu_error, mu_error_d);        // Retiming delay before mult and add
    DelayNUnit #(WIDTH, 1) IN_PIP(clk, reset, x_n_error, x_n_error_d);        // Retiming delay before add

    //assign x_n_error_full = $signed(x_n) * $signed(mu_error_d);
    assign x_n_error_full = $signed(x_n) * $signed(mu_error);
    assign x_n_error_rnd = x_n_error_full + (1<<(QP-1));
    assign x_n_error = x_n_error_rnd[QP+:WIDTH];

    assign new_weight = weight + x_n_error_d;

    always @ ( posedge clk )
    if ( reset )
        weight <= (RESET_VAL <<< QP);
    else
        weight <= new_weight;

endmodule