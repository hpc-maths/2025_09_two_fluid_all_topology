from paraview.simple import *
import os
import subprocess
import xml.etree.ElementTree as ET

overwrite = ("--overwrite" in sys.argv)

xdmf_files = sorted(f for f in os.listdir('.') if f.endswith('.xdmf'))
N_frames = len(xdmf_files)
print(f"Number of .xdmf files: {N_frames}")
if (N_frames==0):
    print("No .xdmf files => run the code and put the files here in order to run this script")
    quit()

pvsm_file = "visu_p2.pvsm"

tree = ET.parse(pvsm_file)
root = tree.getroot()

for prop in root.findall(".//Property[@name='FileNames']"):
    for elem in list(prop.findall("Element")):
        prop.remove(elem)

    for i, filename in enumerate(xdmf_files):
        ET.SubElement(prop, "Element", index=str(i), value=filename)

    prop.set("number_of_elements", str(N_frames))

for prop in root.findall(".//Property[@name='TimestepValues']"):
    # Supprimer les anciens <Element>
    for elem in list(prop.findall("Element")):
        prop.remove(elem)
    # Ajouter les nouveaux (valeurs 0,1,2,...N_frames-1)
    for i in range(N_frames):
        ET.SubElement(prop, "Element", index=str(i), value=str(i))
    # Mettre à jour le nombre d’éléments
    prop.set("number_of_elements", str(N_frames))

# Sauvegarder le nouveau pvsm
tree.write(pvsm_file)
print(f"Updated pvsm with {len(xdmf_files)} .xdmf files")


LoadState(pvsm_file)

renderView = GetActiveViewOrCreate('RenderView')


timeKeeper = GetTimeKeeper()
timeSteps = timeKeeper.TimestepValues


for i, t in enumerate(timeSteps, start=0):
    img_name = f"export_{i:04d}.png"
    if not overwrite and os.path.exists(img_name):
        print(f"Skipping {img_name} (already exists)")
        continue
    print("generating frame " + str(i) + "/"+str(N_frames-1))
    renderView.ViewTime = t
    SaveScreenshot(img_name, renderView)


print("Done generating images")

print("Running ffmpeg to create video...")
ffmpeg_cmd = [
    "ffmpeg",
    "-y",                    
    "-framerate", "5",      
    "-i", "export_%04d.png",
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "2D_Riemann_p2.mp4"
]

print("Running ffmpeg to create video...")
subprocess.run(ffmpeg_cmd, check=True)
print("Video saved as output.mp4")