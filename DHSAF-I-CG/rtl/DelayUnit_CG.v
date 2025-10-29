`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Author: Pavan Kumar
// Create Date: 01-08-2025
// Module name: DelayUnit_CG.v
//////////////////////////////////////////////////////////////////////////////////

module DelayUnit_CG#( parameter LENGTH = 16, WIDTH = 8 )
    ( input clk,
      input reset,
      input [LENGTH-1:0] cg_en,
      input [LENGTH*WIDTH-1:0] reg_in,
      output [LENGTH*WIDTH-1:0] reg_out
    );
    
    wire [WIDTH-1:0] reg_in_unpacked[LENGTH-1:0], reg_out_unpacked[LENGTH-1:0];
    wire clk_gated[LENGTH-1:0];

    // genvar ind;
    // generate
    //     for ( ind = 0; ind < LENGTH; ind=ind+1 )
    //     begin:reg_pack_unpack
            
    //         assign reg_in_unpacked[ind] = reg_in[WIDTH*ind+:WIDTH];
    //         assign clk_gated[ind]       = clk & cg_en[ind];
    //         // CLKAND2X2 cg_and(clk_gated[ind], clk, cg_en[ind]);
            
    //         DelayNUnit#(WIDTH, 1) CG_REG(clk_gated[ind], reset, reg_in_unpacked[ind], reg_out_unpacked[ind]);

    //         assign reg_out[WIDTH*ind+:WIDTH] = reg_out_unpacked[ind];
    //     end
    // endgenerate

    genvar ind;
    generate
        for ( ind = 0; ind < LENGTH; ind=ind+1 )
        begin:reg_pack_unpack
            
            assign reg_in_unpacked[ind] = reg_in[WIDTH*ind+:WIDTH];
            // assign clk_gated[ind]       = clk & cg_en[ind];
            // CLKAND2X2 cg_and(clk_gated[ind], clk, cg_en[ind]);
            
            DelayNUnit_en#(WIDTH, 1) CG_REG(clk, reset, cg_en[ind], reg_in_unpacked[ind], reg_out_unpacked[ind]);

            assign reg_out[WIDTH*ind+:WIDTH] = reg_out_unpacked[ind];
        end
    endgenerate
        
    endmodule
