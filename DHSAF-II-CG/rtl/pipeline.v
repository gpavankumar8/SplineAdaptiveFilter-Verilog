`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 10-08-2022
// Module name: pipeline.v
//////////////////////////////////////////////////////////////////////////////////

module pipeline#( parameter BITSIZE = 8, N = 16 )
    ( input clk,
      input reset,
      input [BITSIZE-1:0] reg_in,
      output [N*BITSIZE-1:0] reg_out_packed
    );
    
    reg [BITSIZE-1:0] shift_reg[N-1:0];
    
    assign reg_out = shift_reg[N-1];
    
    always @ (posedge clk)
    if( reset )
    begin
        shift_reg[0] <= 0;
    end
    else
    begin
        shift_reg[0] <= reg_in;            
    end

    genvar i;
    generate
        for( i = 1; i < N; i = i+1 )
        begin:stage
            always @ (posedge clk)
            if( reset )
            begin
                shift_reg[i] <= 0;            
            end
            else
            begin
                shift_reg[i] <= shift_reg[i-1];            
            end
       end
    endgenerate
        
    genvar ind;
    generate
      
      for ( ind = 0; ind < N; ind=ind+1 )
        begin:outreg
          assign reg_out_packed[BITSIZE*ind+:BITSIZE] = shift_reg[ind];
        end
      
    endgenerate

    endmodule