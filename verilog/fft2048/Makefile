# fft2048 -- Icarus Verilog flow
#   make          generate twiddles, simulate, bit-true check
#   make sim      simulate only
#   make wave     simulate with VCD dump, open GTKWave
#   make model    Python model vs numpy (needs numpy)
#   make clean

IVERILOG ?= iverilog
VVP      ?= vvp
PYTHON   ?= python3

RTL = fft_sdp_ram.v fft_buf_bank.v fft2048_twiddle_rom.v fft2048_bfly.v fft2048.v
TB  = tb_fft2048.v

all: check

tw2048_re.hex tw2048_im.hex: gen_twiddle2048.py
	$(PYTHON) gen_twiddle2048.py

tb_fft2048.vvp: $(RTL) $(TB) tw2048_re.hex tw2048_im.hex
	$(IVERILOG) -g2005 -Wall -o $@ $(TB) $(RTL)

sim: tb_fft2048.vvp
	$(VVP) -n tb_fft2048.vvp

check: sim
	$(PYTHON) fft2048_model.py check tb_vec_*.txt

wave: tb_fft2048.vvp
	$(VVP) -n tb_fft2048.vvp +vcd
	gtkwave tb_fft2048.vcd &

model:
	$(PYTHON) fft2048_model.py selftest

clean:
	rm -f tb_fft2048.vvp tb_fft2048.vcd tb_vec_*.txt

.PHONY: all sim check wave model clean
