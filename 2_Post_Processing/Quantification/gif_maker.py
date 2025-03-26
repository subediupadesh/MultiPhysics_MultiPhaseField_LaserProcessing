import cv2
import imageio
import os
from natsort import natsorted

def create_gif_from_pngs(IP_folder, OP_folder, FPS=30, resize_factor=0.5):
    images = []
    png_files = natsorted([f for f in os.listdir(IP_folder) if f.endswith('.png')])

    for file_name in png_files:
        img_path = os.path.join(IP_folder, file_name)
        img = cv2.imread(img_path)

        if resize_factor != 1.0:
            img = cv2.resize(img, None, fx=resize_factor, fy=resize_factor, interpolation=cv2.INTER_AREA)
        
        # Convert the image from BGR (OpenCV format) to RGB (imageio format)
        img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
        images.append(img_rgb)

    imageio.mimsave(OP_folder, images, format='GIF', fps=FPS,  loop=0)
    print(f"GIF created successfully: {OP_folder}")

m = int(input("Enter \n1 for Model 1, \n2 for Model 2 \n3 for Model 3 \n4 for Model 4 \n:"))
i = int(input("Enter \n1 for Area video, \n2 for Phase Growth video, \n3 for Temperature Distribution video, \n4 for Lead Lag Dynamics video \n:"))


if m ==1:
    model='m1'
elif m==2:
    model='m2'
elif m==3:
    model='m3'
else:
    model='m4'

if i ==1:
    name = 'area'
    animation_name = 'Grain_Growth_(Area).gif'
elif i==2:
    name = 'phases'
    animation_name = 'Phase_Growth.gif'
elif i==3:
    name = 'temperature'
    animation_name = 'Temperature_Distribution.gif'
elif i==4:
    name = 'isotherm'
    animation_name = 'Isotherm.gif'


folder_in = f'animation/{model}/temporary_figures/{name}/'
folder_out = f'animation/gif/{model}/{animation_name}'


create_gif_from_pngs(IP_folder=folder_in, OP_folder=folder_out, FPS=30, resize_factor=0.5)


# area_in = f'/animation/{model}/temporary_figures/{name}/'
# phase_in = f'/animation/{model}/temporary_figures/{name}/'
# isotherm_in = f'/animation/{model}/temporary_figures/{name}/'
# temperature_in = f'/animation/{model}/temporary_figures/{name}/'

# area_out = f'animation/gif/{model}/{animation_name}'
# phase_out = f'animation/gif/{model}/{animation_name}'
# isotherm_out = f'animation/gif/{model}/{animation_name}'
# temperature_out = f'animation/gif/{model}/{animation_name}'

# phase_out = f'../video_animations/{heat_source}/{free_energy}/gif/Phase.gif'
# temp_out = f'../video_animations/{heat_source}/{free_energy}/gif/Temperature.gif'
# vel_out = f'../video_animations/{heat_source}/{free_energy}/gif/Velocity.gif'


# create_gif_from_pngs(IP_folder=temp_in, OP_folder=temp_out, FPS=30, resize_factor=0.5)
# create_gif_from_pngs(IP_folder=vel_in, OP_folder=vel_out, FPS=30, resize_factor=0.5)