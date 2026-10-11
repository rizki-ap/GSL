# xspec_phat -- Icarus Verilog flow
#   make          generate ROM, simulate, bit-true + accuracy check
#   make sim      simulate only
#   make wave     simulate with waveform dump (tb_xspec_phat.fst, open in GTKWave)
#   make clean

IVERILOG ?= iverilog
VVP      ?= vvp
PYTHON   ?= python3
FFT      := ../fft2048

RTL = $(FFT)/fft_sdp_ram.v $(FFT)/fft_buf_bank.v $(FFT)/fft2048_twiddle_rom.v \
      $(FFT)/fft2048_bfly.v $(FFT)/fft2048.v xspec_phat_unit.v xspec_phat.v
TB  = tb_xspec_phat.v

all: check

rsqrt1024.hex: gen_rsqrt_rom.py
	$(PYTHON) gen_rsqrt_rom.py

$(FFT)/tw2048_re.hex $(FFT)/tw2048_im.hex:
	cd $(FFT) && $(PYTHON) gen_twiddle2048.py

tb_xspec_phat.vvp: $(RTL) $(TB) rsqrt1024.hex $(FFT)/tw2048_re.hex
	$(IVERILOG) -g2005 -Wall -o $@ $(TB) $(RTL)

sim: tb_xspec_phat.vvp
	$(VVP) -n tb_xspec_phat.vvp

check: sim
	$(PYTHON) xspec_phat_model.py check tb_xp_*.txt

wave: tb_xspec_phat.vvp
	$(VVP) -n tb_xspec_phat.vvp +vcd +fst -fst
	@echo "Waveform: tb_xspec_phat.fst  (open with: gtkwave tb_xspec_phat.fst)"

clean:
	rm -f tb_xspec_phat.vvp tb_xspec_phat.vcd tb_xspec_phat.fst tb_xp_*.txt

.PHONY: all sim check wave clean
