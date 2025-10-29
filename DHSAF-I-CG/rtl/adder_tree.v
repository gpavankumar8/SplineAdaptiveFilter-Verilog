`timescale 1ns / 1ps

//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 11-08-2022
// Module name: adder_tree.v
//
// Modified by: Vishnu P S
// Stage count starts from 0. Eg: for LEN=1024 -> Number of stages needed = 10; Stages will be from 0 till 9.
// Best stage at which pipeline can be done is found out experimentally.
// Accordingly, it is found that for WIDTH = 16, a pipeline stage has to be added at every 4th stage.
// Therefore a Pipeline stage will be added at every {piping_start_stage+1}th stage, starting from {piping_start_stage}th stage
//////////////////////////////////////////////////////////////////////////////////

module adder_tree
#(parameter WIDTH=16, LEN=1024, piping_start_stage=3)
(
    input [LEN*WIDTH-1:0] adder_tree_in_packed,
    input clk,reset,
    output [WIDTH-1:0] adder_tree_out
);
    wire signed [WIDTH-1:0] stage[$clog2(LEN)-1:0][LEN/2-1:0];
    wire signed [WIDTH-1:0] adder_tree_in[LEN-1:0];
    parameter max_possible_stage = $clog2(LEN)-1;
    
    // Unpack adder_tree_in
    genvar i;
    generate
        for ( i = 0; i < LEN; i=i+1 )
        begin:unpack
            assign adder_tree_in[i] = adder_tree_in_packed[WIDTH*i+:WIDTH];
        end
    endgenerate

    // Stage 0 of adder tree
    genvar j;
    generate
    for(j = 0; j < LEN/2; j = j+1)
        begin : tree0
            assign stage[0][j] = adder_tree_in[2*j] + adder_tree_in[2*j+1];
        end
    endgenerate 

    // Rest of the adder stages
    genvar p,q,r,s,t,u,v,w,x;
    generate
        // Checks whether pipeline is needed or not
        if (max_possible_stage-piping_start_stage > 0) 
        begin    

            // Subsequent stages until first pipeline stage
            for (p = 1; p <= piping_start_stage; p = p+1)
            begin : tree_before_pipeline_stage
                for(q = 0; q < LEN/4; q = q+1)
                begin : tree_before_pipeline_stage_inner
                    assign stage[p][q] = stage[p-1][2*q] + stage[p-1][(2*q)+1];
                end
            end

            // Stages from first pipeline till last
            for(r = piping_start_stage; r < max_possible_stage; r = r+(piping_start_stage+1)) 
            begin: multipipe

                // Pipeline register array
                wire [WIDTH-1:0] pip_reg[(LEN/2**(r+1))-1:0]; 

                // Pushing into pipeline reg
                for (s=0; s < ((LEN/2**(r+1))); s=s+1)
                begin: pipe_in 
                    //register #(.data_width(WIDTH)) piping(.inData(stage[r][s]),.outData(pip_reg[s]),.clk(clk),.reset(reset));
                    DelayNUnit #(WIDTH, 1) piping(clk, reset, stage[r][s], pip_reg[s]);
                end

                // The stage immediately after pipelining
                for (t=0; t < (LEN/2**(r+1))/2; t=t+1)
                begin: pipe_out 
                    assign stage[r+1][t] = pip_reg[2*t] + pip_reg[(2*t)+1];
                end
            
                // Subsequent stages after pipelined stage until next pipeline stage (or) last stage whichever comes first
                for(u = r+2; (u <= max_possible_stage) && (u <= r+(piping_start_stage+1)); u = u+1)
                begin : tree_after_pipeline_stage
                    for(v = 0; v < LEN/4; v = v+1)
                    begin : tree_after_pipeline_stage_inner
                        assign stage[u][v] = stage[u-1][2*v] + stage[u-1][(2*v)+1];
                    end
                end

            end
        end

        else 
        begin
        // This section works if no pipelining is needed

            // Subsequent stages after stage 0
            for(w = 1; (w <= max_possible_stage); w = w+1)
            begin : tree_with_no_pipeline
                for(x = 0; x < LEN/4; x = x+1)
                begin : tree_with_no_pipeline_inner
                    assign stage[w][x] = stage[w-1][2*x] + stage[w-1][2*x+1];
                end
            end

        end
    endgenerate 

    assign adder_tree_out = stage[max_possible_stage][0];

endmodule