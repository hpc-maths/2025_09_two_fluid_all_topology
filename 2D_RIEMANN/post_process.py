# save this as render_xdmf.py
from paraview.simple import *
import glob
import os
import subprocess
# Folder containing your XDMF files
xdmf_folder = "RESULTS"
output_folder = "frames"
# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)
# List all XDMF files sorted by time
xdmf_files = sorted(glob.glob(os.path.join(xdmf_folder, "*.xdmf")))
print(xdmf_files)
# Loop over each XDMF file
for i, xdmf_file in enumerate(xdmf_files):
    if(i>20):
        break
    print(f"Processing {xdmf_file}")
    # Load the XDMF file
    reader = XDMFReader(FileNames=[xdmf_file])
    reader.UpdatePipeline() 
    
    renderView = GetActiveViewOrCreate('RenderView')
    display = Show(reader, renderView)
    ColorBy(display, ('CELLS', 'rho1'))
    ShowScalarBar = True
    display.SetScalarBarVisibility(renderView, ShowScalarBar)
    Render()
    scalarBar = GetScalarBar(display, renderView)
    if scalarBar is not None:
        scalarBar.Title = "density 1"
        scalarBar.ComponentTitle = ""
        scalarBar.Visibility = 1
    # Reset camera to fit the data
    renderView.ResetCamera()
    # Save screenshot
    frame_filename = os.path.join(output_folder, f"frame_{i:04d}.png")
    SaveScreenshot(frame_filename, renderView, ImageResolution=[1280, 720])
    
    # Clean up for next iteration
    Delete(reader)
    del reader
print("All frames saved.")


ffmpeg_cmd = [
    "ffmpeg",
    "-y",                    # overwrite output if exists
    "-framerate", "10",      # change 10 → your desired FPS
    "-i", "frames/frame_%04d.png",
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "2D_Riemann_rho_1.mp4"
]

print("Running ffmpeg to create video...")
subprocess.run(ffmpeg_cmd, check=True)
print("Video saved as output.mp4")