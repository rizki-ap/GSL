# gcc_engine -- Milestone A (Icarus Verilog flow)
#   make            prepare jobs -> simulate -> check  (full Milestone A)
#   make prepare    build jobs (synthetic + your generator in ../../python)
#   make sim        simulate only (needs prepared jobs)
#   make check      Python check of the last simulation
#   make wave       simulate with waveform dump (tb_gcc_engine.fst)
#   make clean      remove generated files (keeps nothing from prepare)

IVERILOG ?= iverilog
VVP      ?= vvp
PYTHON   ?= python3
FFT      := ../fft2048
XSP      := ../xspec_phat

RTL = $(FFT)/fft_sdp_ram.v $(FFT)/fft_buf_bank.v $(FFT)/fft2048_twiddle_rom.v \
      $(FFT)/fft2048_bfly.v $(FFT)/fft2048.v \
      $(XSP)/xspec_phat_unit.v $(XSP)/xspec_phat.v \
      cap_ring.v gcc_loader.v peak_search.v gcc_engine.v
TB  = tb_gcc_engine.v

all: prepare sim check

prepare:
	$(PYTHON) milestone_a.py prepare

tb_gcc_engine.vvp: $(RTL) $(TB)
	$(IVERILOG) -g2005 -Wall -o $@ $(TB) $(RTL)

sim: tb_gcc_engine.vvp
	$(VVP) -n tb_gcc_engine.vvp

check:
	$(PYTHON) milestone_a.py check

wave: tb_gcc_engine.vvp
	$(VVP) -n tb_gcc_engine.vvp +vcd +fst -fst
	@echo "Waveform: tb_gcc_engine.fst  (open with: gtkwave tb_gcc_engine.fst)"

clean:
	rm -rf tb_gcc_engine.vvp tb_gcc_engine.vcd tb_gcc_engine.fst \
	       ma_jobs.txt ma_ring_*.hex ma_prep.json ma_out_*.txt gen_work

.PHONY: all prepare sim check wave clean
