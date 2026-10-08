# ModelSim / Questa (Intel FPGA edition)
#   cd verilog/fft2048
#   python3 gen_twiddle2048.py
#   vsim -c -do do_fft2048.do
vlib work
vlog fft_sdp_ram.v fft_buf_bank.v fft2048_twiddle_rom.v fft2048_bfly.v fft2048.v tb_fft2048.v
vsim -voptargs=+acc work.tb_fft2048
add wave -noupdate /tb_fft2048/clk
add wave -noupdate /tb_fft2048/dut/state
add wave -noupdate -radix unsigned /tb_fft2048/dut/stage
add wave -noupdate -radix unsigned /tb_fft2048/dut/ctr
add wave -noupdate -radix unsigned /tb_fft2048/dut/sh_cur
add wave -noupdate -radix unsigned /tb_fft2048/dut/exp_acc
add wave -noupdate /tb_fft2048/dut/fft_done
add wave -noupdate -radix unsigned /tb_fft2048/dut/fft_bfp_exp
add wave -noupdate /tb_fft2048/dut/fft_err
run -all
