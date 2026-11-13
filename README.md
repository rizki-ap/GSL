# Gunshot-Locator
research about implementation of Gunshot Locator 

Product Reference : 
1. PILAR-V. Link :
2. Boomerang

Publication Reference : 
- main : https://pub.dega-akustik.de/DAGA_1999-2008/data/articles/001903.pdf
- other in /reference

Gunshot Recording Dataset :
- https://zenodo.org/record/7004819 ; Publication : https://doi.org/10.1016/j.dib.2023.109091
- https://cadreforensics.com/audio/
- data in /sound (Cipatat, Cipatat_clean, NIJ, Zenodo, Cadre; see sound/readme.md)

Target Implementation Component Data :
- in /devices

Python Algorithm : (/python)
## Generator Part
- gen_config.ini : all parameter needed to configure generation of real gunshot signal
- gs_gen_physic.py : physics parameter setup, 
- gs_gen_clean_signal.py : generate clean gunshot signal (shockwave, muzzle blast, optional ground reflection) in multichannel wav file and *.json 
- gs_gen_add_noise.py : adding environment (noise layers, echoes, reverb), create multichannel wav file and *.json
- gs_gen_apply_adc.py : adding noise in mics, afe and adc.

## Detector Part
- det_config.ini : all parameter needed to configure detection gunshot shooter origin
- gs_det_signal_prepare.py :detecting sw & mb peak, normalization
- gs_det_tdoa.py : cross-correlation (default) or gcc-phat, create 6 pair tDoA
- gs_det_shooter_locator.py : bearing from tDoA, drag-aware range, determine shooter location
- gs_det_classify_bullet.py : determine bullet property, mach, caliber, etc

## Utility 
- util_get_accoustic_param.py : extracting accoustic parameter, include noise from wav file
- util_plot_signal.py : plotting signal in multichannel wav file
- util_signal_similarity.py : comparing measured and synthesized gunshots
- util_gen_cipatat_clean.py : clean synthetic twins of the Cipatat recordings
- full list and details: python/project.md

FPGA rtl implementation :
- in /verilog

