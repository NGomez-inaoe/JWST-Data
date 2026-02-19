from pathlib import Path
import shutil


#this moves subdirectories and files from path to final_path
"""
path=Path("./MAST/JWST/")
final_path=Path("./MAST/")

for item in path.iterdir():
    if item.is_dir():
        shutil.move(str(item), final_path / item.name)
"""
     
#this moves all the 1d files from their folders to the general folder for MAST
source=Path("./MAST")
destination=Path("./Spectra1D/MAST")

destination.mkdir(parents=True, exist_ok=True)

for subdir in source.iterdir():
    for file in subdir.glob("*x1d.fits"):
        shutil.move(str(file), destination / file.name)
        
        
source_jades=Path("./JADES")
destination_jades_1d=Path("./Spectra1D/JADES")

destination_jades_1d.mkdir(parents=True, exist_ok=True)

for subdir in source_jades.iterdir():
    for file in subdir.glob("*x1d.fits"):
        shutil.move(str(file), destination_jades_1d / file.name)
        
#same for 2d files
source_mast=Path("./MAST")
destination_mast_2d=Path("./Spectra2D/MAST")

destination_mast_2d.mkdir(parents=True, exist_ok=True)

for subdir in source_mast.iterdir():
    for file in subdir.glob("*s2d.fits"):
        shutil.move(str(file), destination_mast_2d / file.name)
        
        
source_jades=Path("./JADES")
destination_jades_2d=Path("./Spectra2D/JADES")

destination_jades_2d.mkdir(parents=True, exist_ok=True)

for subdir in source_jades.iterdir():
    for file in subdir.glob("*s2d.fits"):
        shutil.move(str(file), destination_jades_2d / file.name)

