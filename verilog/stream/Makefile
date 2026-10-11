# stream -- Milestone B (continuous detection chain)
#   make              prepare scenarios -> build -> simulate all -> check
#   make prepare      coefficients + scenarios from your generator (../../python)
#   make sim          simulate every scenario in mb_list.txt
#   make check        Python bit-true + accuracy check
#   make SIM=icarus   use Icarus instead of Verilator (much slower; FRAME_DIV 100)
#   make FRAME_DIV=1000  real-time cadence (default 250: 4x faster, identical
#                     results, and a harsher test of the scheduler)
#   make clean

SIM      ?= verilator
FRAME_DIV ?= 250
PYTHON   ?= python3
FFT      := ../fft2048
XSP      := ../xspec_phat
ENG      := ../gcc_engine

RTL = $(FFT)/fft_sdp_ram.v $(FFT)/fft_buf_bank.v $(FFT)/fft2048_twiddle_rom.v \
      $(FFT)/fft2048_bfly.v $(FFT)/fft2048.v \
      $(XSP)/xspec_phat_unit.v $(XSP)/xspec_phat.v \
      $(ENG)/cap_ring.v $(ENG)/gcc_loader.v $(ENG)/peak_search.v $(ENG)/gcc_engine.v \
      replay_src.v timebase.v fe_cond.v mb_decim.v sw_det.v mb_det.v evt_sched.v gsl_core.v
TB  = tb_gsl_core.v

ifeq ($(SIM),icarus)
  SIMBIN := tb_gsl_core.vvp
  RUN    := vvp -n $(SIMBIN)
else
  SIMBIN := obj_dir/Vtb_gsl_core
  RUN    := $(SIMBIN)
endif

all: prepare sim check

prepare:
	$(PYTHON) milestone_b.py prepare

obj_dir/Vtb_gsl_core: $(RTL) $(TB)
	verilator --binary --timing -O3 -Wno-fatal -Wno-lint -Wno-style -GFRAME_DIV=$(FRAME_DIV) \
	    --top-module tb_gsl_core -o Vtb_gsl_core $(TB) $(RTL)

tb_gsl_core.vvp: $(RTL) $(TB)
	iverilog -g2005 -Ptb_gsl_core.FRAME_DIV=100 -o $@ $(TB) $(RTL)

sim: $(SIMBIN)
	@while read n f; do \
	  echo "--- $$n ---"; \
	  $(RUN) +stim=mb_$$n.hex +nframes=$$f +cfg=mb_cfg.txt +out=mb_$$n | grep -v "frames replayed"; \
	done < mb_list.txt

check:
	$(PYTHON) milestone_b.py check

clean:
	rm -rf obj_dir tb_gsl_core.vvp mb_*.hex mb_*_truth.json mb_*_x.txt mb_*_y.txt \
	       mb_*_ev.txt mb_cfg.txt mb_list.txt gen_work __pycache__ build.log sim.log

.PHONY: all prepare sim check clean
