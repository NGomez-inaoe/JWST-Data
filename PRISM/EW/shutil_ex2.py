from pathlib import Path
import shutil

# Carpetas originales
scripts_dir = Path("EW/EW-Scripts")
plots_dir = Path("EW/EW-Plots")
data_dir = Path("EW/EW-Data-Output")

# Carpeta destino
objects_dir = Path("EW/EW-Objects")
objects_dir.mkdir(exist_ok=True)

# Iterar sobre los scripts
for script_file in scripts_dir.glob("Equivalent_Widths_*.py"):
    
    # Extraer xxx
    name = script_file.stem.replace("Equivalent_Widths_", "")
    
    # Crear subcarpeta EW-Objects/xxx
    target_subdir = objects_dir / name
    target_subdir.mkdir(exist_ok=True)
    
    # Construir nombres correspondientes
    pdf_file = plots_dir / f"EW-{name}.pdf"
    tsv_file = data_dir / f"EW_output_{name}.tsv"
    
    # Copiar script
    shutil.copy2(script_file, target_subdir / script_file.name)
    
    # Copiar PDF si existe
    if pdf_file.exists():
        shutil.copy2(pdf_file, target_subdir / pdf_file.name)
    else:
        print(f"No se encontró {pdf_file.name}")
    
    # Copiar TSV si existe
    if tsv_file.exists():
        shutil.copy2(tsv_file, target_subdir / tsv_file.name)
    else:
        print(f"No se encontró {tsv_file.name}")

print("Estructura creada correctamente.")