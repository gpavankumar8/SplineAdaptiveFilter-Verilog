`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 20.01.2021 12:25:06
// Design Name: 
// Module Name: hsaf_tb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


`define clock 10

module hsaf_tb();

  parameter L_ORD = 8;
  parameter Q_ORD = 4;
  parameter WIDTH = 16;
  parameter QP = 12;
  parameter N = 20000;
  parameter NUM_TRIAL = 50;
  parameter DelX_inv = 2;
  parameter Q = 13;
  parameter RET = 4;

  reg              clk, rst;
  reg  [WIDTH-1:0] x, d;
  wire [WIDTH-1:0] y, error;

  hsaflms_top#(L_ORD, Q_ORD, WIDTH, QP, DelX_inv, Q) HSAF_DUT( .clk(clk), .reset(rst), .signal_in(x), .desired_in(d), .filter_out_d(y), .error_d(error));
  // hsaflms_top HSAF_DUT( .clk(clk), .reset(rst), .signal_in(x), .desired_in(d), .filter_out_d(y), .error_d(error));

  //initial begin
    //$dumpfile("dump.vcd");
    //$dumpvars(0);
  //end
  
  // Clock generation
  initial begin
    clk=1'b1;
    forever #(`clock/2) clk=~clk; //Clock Generator
  end
  
  // Test bench related variables
  reg [999:0]     fname_x, fname_d, fname_err_rtl;
  reg [WIDTH-1:0] x_val[0:N-1], d_val[0:N-1];
  integer i, k, f, trial, ret;

  initial begin

    for(trial = 1; trial <= NUM_TRIAL; trial = trial + 1)
    begin

        // File name generation and setup
        $sformat(fname_x, "./inputs/L8_v5/x_values%0d.txt",trial);
        $sformat(fname_d, "./inputs/L8_v5/d_values%0d.txt",trial);
        $sformat(fname_err_rtl, "./outputs/L8_v5/error_rtl%0d.txt",trial);
        f = $fopen(fname_err_rtl,"w");

        // Read x and d values from file
        $readmemh(fname_x, x_val);
        $readmemh(fname_d, d_val);

        ret = RET;
        
        // Driving stimulus
        rst = 1'b1;
        
        #10 d = d_val[0];
            x = x_val[0];
        #20 rst = 1'b0;   
        
        for ( i=1; i<N; i=i+1)
        begin      
            @(posedge clk)
            begin
              #5
              if (ret == 0)
                $fwrite(f,"%b\n",error);
              else
                ret = ret - 1;

              d = d_val[i];
              x = x_val[i];
            end

        end

        for ( i=0; i<RET; i=i+1)
        begin      
          @(posedge clk)
          begin
              #5 $fwrite(f,"%b\n",error);
          end
        end

        #10  $fwrite(f,"%b\n",error);
        #10  $fwrite(f,"%b\n",error);
             $fclose(f);
    end

    #100 rst = 1'b1;
    #20  $finish;
  end

endmodule
