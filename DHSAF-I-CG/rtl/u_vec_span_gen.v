`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 09-05-2023
// Module name: u_vec_span_gen.v
//////////////////////////////////////////////////////////////////////////////////

module u_vec_span_gen
#(parameter WIDTH = 16, Q_ORD = 4, QP = 12, Q = 13, DelX_inv = 2, SPAN_WIDTH = 5)
(
    input clk,
    input reset,
    input [WIDTH-1:0] x_in,
    output [Q_ORD*WIDTH-1:0] u_vec_C_packed,
    output [SPAN_WIDTH-1:0] span_ind
);

    
    wire signed [WIDTH-1:0] u, u_d, u_2d, u2, u2_d, u2_2d, u3, u3_d, x_div_DelX, u_vec_C[Q_ORD-1:0];
    wire signed [2*WIDTH-1:0] u2_full, u3_full, u2_rnd, u3_rnd;

    // Generate u and span index (j)
    assign x_div_DelX = x_in <<< DelX_inv;
    assign u = x_div_DelX[QP-1:0];
    assign span_ind = $signed(x_div_DelX[WIDTH-1:QP]) + $signed((Q-1)/2);
    DelayNUnit #(WIDTH, 1) U_PIP(clk, reset, u, u_d);
    // DelayNUnit #(WIDTH, 1) U_D_PIP(clk, reset, u_d, u_2d);

    // Generate u_vec_C vector
    assign u2_full = u  * u;
    assign u2_rnd = u2_full + (1'b1 << (QP-1));
    assign u2 = u2_rnd[QP+:WIDTH];
    DelayNUnit #(WIDTH, 1) U2_PIP(clk, reset, u2, u2_d);
    // DelayNUnit #(WIDTH, 1) U2_D_PIP(clk, reset, u2_d, u2_2d);

    assign u3_full = u2_d * u_d;
    assign u3_rnd = u3_full + (1'b1 << (QP-1));
    assign u3 = u3_rnd[QP+:WIDTH];
    // DelayNUnit #(WIDTH, 1) U3_PIP(clk, reset, u3, u3_d);

    // Hardcoded u_vec_C for CR spline matrix and P = 3 (Q_ORD = 4)
    // assign u_vec_C[0] = (-1*u3 + 2*u2 - u) >>> 1;
    // assign u_vec_C[1] = ( 3*u3 - 5*u2 + (2'b10 << QP)) >>> 1;
    // assign u_vec_C[2] = (-3*u3 + 4*u2 + u) >>> 1;
    // assign u_vec_C[3] = (   u3 -   u2    ) >>> 1;

    assign u_vec_C[0] = (-1*u3 + 2*u2_d - u_d) >>> 1;
    assign u_vec_C[1] = ( 3*u3 - 5*u2_d + (2'b10 << QP)) >>> 1;
    assign u_vec_C[2] = (-3*u3 + 4*u2_d + u_d) >>> 1;
    assign u_vec_C[3] = (   u3 -   u2_d       ) >>> 1;

    // Packing u_vec_C
    genvar ind;
    generate
        for ( ind = 0; ind < Q_ORD; ind=ind+1 )
        begin:u_vec_C_pack
            assign u_vec_C_packed[WIDTH*ind+:WIDTH] = u_vec_C[ind];
        end
    endgenerate    

endmodule