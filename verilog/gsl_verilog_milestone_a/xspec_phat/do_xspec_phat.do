# ModelSim / Questa (Intel FPGA edition)
#   cd verilog/xspec_phat
#   vsim -c -do do_xspec_phat.do
#   python3 xspec_phat_model.py check tb_xp_*.txt
vlib work
vlog ../fft2048/fft_sdp_ram.v ../fft2048/fft_buf_bank.v ../fft2048/fft2048_twiddle_rom.v \
     ../fft2048/fft2048_bfly.v ../fft2048/fft2048.v \
     xspec_phat_unit.v xspec_phat.v tb_xspec_phat.v
vsim -voptargs=+acc work.tb_xspec_phat
add wave -noupdate /tb_xspec_phat/clk
add wave -noupdate /tb_xspec_phat/u_xp/state
add wave -noupdate -radix unsigned /tb_xspec_phat/u_xp/slot
add wave -noupdate -radix unsigned /tb_xspec_phat/u_xp/ph
add wave -noupdate /tb_xspec_phat/u_xp/w_we
add wave -noupdate -radix unsigned /tb_xspec_phat/u_xp/w_waddr
add wave -noupdate /tb_xspec_phat/u_fft/fft_done
add wave -noupdate -radix unsigned /tb_xspec_phat/u_fft/fft_bfp_exp
run -all
