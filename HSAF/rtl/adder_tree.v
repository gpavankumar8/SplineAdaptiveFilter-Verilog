`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 11-08-2022
// Module name: adder_tree.v
//////////////////////////////////////////////////////////////////////////////////

module adder_tree
#(parameter WIDTH=16, LEN=64)
(
    input [LEN*WIDTH-1:0] adder_tree_in_packed,
    output [WIDTH-1:0] adder_tree_out
);

    wire signed [WIDTH-1:0] stage[$clog2(LEN)-1:0][LEN/2-1:0];
    wire signed [WIDTH-1:0] adder_tree_in[LEN-1:0];
    
    // Unpack adder_tree_in
    genvar ind;
    generate
        for ( ind = 0; ind < LEN; ind=ind+1 )
        begin:outreg
            assign adder_tree_in[ind] = adder_tree_in_packed[WIDTH*ind+:WIDTH];
        end

    endgenerate

    // Stage 0 adder tree
    genvar k;
    generate
    for(k = 0; k < LEN/2; k = k+1)
        begin : tree0
            assign stage[0][k] = adder_tree_in[2*k] + adder_tree_in[2*k+1];
        end
    endgenerate 

    // Subsequent stages using genvar
    genvar i,j;
    generate
    for(i = 1; i < $clog2(LEN); i = i+1)
    begin : tree
          
        for(j = 0; j < LEN/4; j = j+1)
        begin : tree_inner
            assign stage[i][j] = stage[i-1][2*j] + stage[i-1][2*j+1];
        end
    end
    endgenerate 

    assign adder_tree_out = stage[$clog2(LEN)-1][0];

endmodule