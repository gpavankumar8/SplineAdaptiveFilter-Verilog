`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 10-08-2022
// Module name: DelayNUnit_en.v
//////////////////////////////////////////////////////////////////////////////////

module DelayNUnit_en#( parameter BITSIZE = 8, N = 16  )
    ( input clk,
      input reset,
      input enable,
      input [BITSIZE-1:0] reg_in,
      output [BITSIZE-1:0] reg_out
    );
    
    reg [BITSIZE-1:0] shift_reg[N-1:0];
    
    assign reg_out = shift_reg[N-1];
    
    always @ (posedge clk)
    if( reset )
    begin
        shift_reg[0] <= 0;
    end
    else if( enable )
    begin
        shift_reg[0] <= reg_in;            
    end
    else
    begin
        shift_reg[0] <= shift_reg[0];     
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
            else if( enable )
            begin
                shift_reg[i] <= shift_reg[i-1];            
            end
            else
            begin
                shift_reg[i] <= shift_reg[i];     
            end
        
       end
    endgenerate
        
    endmodule
