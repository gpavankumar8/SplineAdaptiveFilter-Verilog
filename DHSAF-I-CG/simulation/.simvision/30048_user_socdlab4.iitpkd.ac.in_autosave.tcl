
# NC-Sim Command File
# TOOL:	ncsim(64)	15.20-s051
#

set tcl_prompt1 {puts -nonewline "ncsim> "}
set tcl_prompt2 {puts -nonewline "> "}
set vlog_format %h
set vhdl_format %v
set real_precision 6
set display_unit auto
set time_unit module
set heap_garbage_size -200
set heap_garbage_time 0
set assert_report_level note
set assert_stop_level error
set autoscope yes
set assert_1164_warnings yes
set pack_assert_off {}
set severity_pack_assert_off {note warning}
set assert_output_stop_level failed
set tcl_debug_level 0
set relax_path_name 1
set vhdl_vcdmap XX01ZX01X
set intovf_severity_level ERROR
set probe_screen_format 0
set rangecnst_severity_level ERROR
set textio_severity_level ERROR
set vital_timing_checks_on 1
set vlog_code_show_force 0
set assert_count_attempts 1
set tcl_all64 false
set tcl_runerror_exit false
set assert_report_incompletes 0
set show_force 1
set force_reset_by_reinvoke 0
set tcl_relaxed_literal 0
set probe_exclude_patterns {}
set probe_packed_limit 4k
set probe_unpacked_limit 16k
set assert_internal_msg no
set svseed 1
set assert_reporting_mode 0
alias iprof profile
database -open -shm -into waves.shm waves -default
probe -create -database waves hsaf_tb_5_2ns.clk hsaf_tb_5_2ns.d hsaf_tb_5_2ns.error hsaf_tb_5_2ns.f hsaf_tb_5_2ns.fname_d hsaf_tb_5_2ns.fname_err_rtl hsaf_tb_5_2ns.fname_x hsaf_tb_5_2ns.i hsaf_tb_5_2ns.k hsaf_tb_5_2ns.ret hsaf_tb_5_2ns.rst hsaf_tb_5_2ns.trial hsaf_tb_5_2ns.x hsaf_tb_5_2ns.y
probe -create -database waves hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1165.A hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1165.Y hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1167.A hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1167.Y hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1169.A hsaf_tb_5_2ns.HSAF_DUT.EC_sub_16_31_g1169.Y hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][0] }.CK hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][0] }.D hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][0] }.Q hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][0] }.QBINT hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][1] }.CK hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][1] }.D hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][1] }.Q hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][1] }.QBINT hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][2] }.CK hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][2] }.D hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][2] }.Q hsaf_tb_5_2ns.HSAF_DUT.@{\FILT_OUT_PIP_shift_reg_reg[0][2] }.QBINT

simvision -input /home/user/pavan/projects/adaptfilt/ASIC/Retimed/HSAFLMS-lessPIP/simulation/.simvision/30048_user_socdlab4.iitpkd.ac.in_autosave.tcl.svcf
