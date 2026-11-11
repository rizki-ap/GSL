Sound Dataset from multiple gunshot for Algorithm development

- https://cadreforensics.com/audio/download/
- https://zenodo.org/records/7004819
- https://data.jrc.ec.europa.eu/dataset/011799a6-80c1-4416-ae15-17116e2fc32f#description

From Zenodo:


glock_a283_*.wav --- 2021-11-15T11:44:39.555918-05:00; Raspberry Pi	Zero W;	MiniDSP UMA-8 USB Mic Array V2;	Glock 17	9mm;	28.3425767;	-80.78142127; 50;	NORTH	EAST

ruger_04_1s_*.wav --- 2021-11-22T09:59:51.692048-05:00; Raspberry Pi	Zero W;	MiniDSP UMA-8 USB Mic Array V2;	Ruger 556	.223mm;	28.34258;	-80.78142; 50;	NORTH	EAST

ruger_1b_*.wav --- 2021-11-15T12:07:18.185684-05:00; Raspberry Pi	Zero W;	MiniDSP UMA-8 USB Mic Array V2;	Ruger 556	.223mm;	28.3425767;	-80.78142127;	50;	NORTH	EAST



- collection setup : https://www.sciencedirect.com/science/article/pii/S235234092300210X


## Folders (Nov 2026)

| Folder | Contents |
|---|---|
| `Cipatat/` | PILAR-V test-range recordings, 4 ch, 51.2 kHz, 5.56 mm to 120 mm at 88-1470 m. See `Cipatat/metadata.md` |
| `Cipatat_clean/` | Clean (noise- and reverberation-free) synthetic twins of the Cipatat development recordings. See `Cipatat_clean/metadata.md` |
| `Cadre/` | Sample of the Gunshot Audio Forensics Dataset (cadreforensics.com): M16 shots, experiment 17 A/B. See `Cadre/metadata.md` |
| `Zenodo/` | Edge-collected UMA-8 recordings (zenodo.org/records/7004819): 7 channels + mean per shot, rifle, pistol and shotgun. See `Zenodo/metadata.md` |
| `NIJ/` | AR-15 and .308 shots at 500 kHz, 12 mics on a 3 m semicircle around the gun (NIJ grant 2014-DN-BX-K034). See `NIJ/metadata.md` |

