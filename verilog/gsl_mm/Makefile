# gsl_mm -- Milestone C1 (processor interface: registers, mailbox, replay DMA)
#   make            build, then run every Milestone B scenario through gsl_mm and
#                   compare with Milestone B (needs ../stream prepared: make -C ../stream prepare)
#   make build      Verilator build only
#   make clean

FRAME_DIV ?= 250
PYTHON    ?= python3
FFT := ../fft2048
XSP := ../xspec_phat
ENG := ../gcc_engine
STR := ../stream

RTL = $(FFT)/fft_sdp_ram.v $(FFT)/fft_buf_bank.v $(FFT)/fft2048_twiddle_rom.v \
      $(FFT)/fft2048_bfly.v $(FFT)/fft2048.v \
      $(XSP)/xspec_phat_unit.v $(XSP)/xspec_phat.v \
      $(ENG)/cap_ring.v $(ENG)/gcc_loader.v $(ENG)/peak_search.v $(ENG)/gcc_engine.v \
      $(STR)/replay_src.v $(STR)/timebase.v $(STR)/fe_cond.v $(STR)/mb_decim.v \
      $(STR)/sw_det.v $(STR)/mb_det.v $(STR)/evt_sched.v $(STR)/gsl_core.v \
      gsl_mbox.v gsl_replay_dma.v gsl_mm.v
TB  = tb_gsl_mm.v

all: build run

build: obj_dir/Vtb_gsl_mm

obj_dir/Vtb_gsl_mm: $(RTL) $(TB)
	verilator --binary --timing -O3 -Wno-fatal -Wno-lint -Wno-style -GFRAME_DIV=$(FRAME_DIV) \
	    --top-module tb_gsl_mm -o Vtb_gsl_mm $(TB) $(RTL)

run: obj_dir/Vtb_gsl_mm
	$(PYTHON) c1_check.py

clean:
	rm -rf obj_dir c1_*_x.txt c1_*_y.txt c1_*_ev.txt __pycache__ build.log

.PHONY: all build run clean
