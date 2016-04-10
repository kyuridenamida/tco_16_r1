#!C:\Users\kyuridenamida\AppData\Local\Programs\Python\Python35\python.exe
import matplotlib.pyplot as plt
import pygame
from pygame.locals import *
import sys
 
'''
仕様:
引数に指定されている命令列を改行区切りで実行します．
[円]C x y 半径 
[点]P x y
[直線] L x1 y1 x2 y2
[多角形] G x1 y1 x2 y2 x3 y3 ... xn yn
[シーンの終わり] END
'''

scenes = []
current_list = []
for line in sys.stdin:
	line = line.strip()
	if line == "END":
		scenes.append(current_list)
		current_list = []
	else:
		current_list.append(line)

SCREEN_SIZE = (1024, 1024)
 
pygame.init()
screen = pygame.display.set_mode(SCREEN_SIZE)

 
 
pos = 0

while True:
	
	screen.fill((255,255,255)) 
	screen.set_alpha(128)
	pygame.display.set_caption("%d/%d"%(pos+1,len(scenes)))
	for scene in scenes[pos]:
		scene = scene.split(' ')
		type = scene[0]
		args = list(map(int,scene[1:]))
		color = (args[0],args[1],args[2])
		args = args[3:]
		if scene[0] == 'C':
			if args[2] != 0:
				pygame.draw.circle(screen, color, (args[0],args[1]), args[2],1)
		elif scene[0] == 'L':
			vec = (args[2] - args[0], args[3] - args[1])
			p1 = (args[0] + 1024 * vec[0],args[1] + 1024 * vec[1])
			p2 = (args[0] - 1024 * vec[0],args[1] - 1024 * vec[1])
			
			pygame.draw.line(screen, color, p1,p2 )
		elif scene[0] == 'G':
			# poly = pygame.Surface((100,100), pygame.SRCALPHA, 32)
			pygame.draw.polygon(screen, color, list(zip(args[0::2],args[1::2])) )
			# screen.blit(poly, (100,100))
	
	for scene in scenes[pos]:
		scene = scene.split(' ')
		type = scene[0]
		args = list(map(int,scene[1:]))
		color = (args[0],args[1],args[2])
		args = args[3:]
		if scene[0] == 'C':
			if args[2] != 0:
				pygame.draw.circle(screen, color, (args[0],args[1]), args[2],1)
		elif scene[0] == 'G':
			# poly = pygame.Surface((100,100), pygame.SRCALPHA, 32)
			pygame.draw.polygon(screen, color, list(zip(args[0::2],args[1::2])),1 )
			# screen.blit(poly, (100,100))
	
	pygame.display.update()
	
	for event in pygame.event.get():
		if event.type == QUIT:
			sys.exit()
		elif event.type == KEYUP:
			if event.key == K_LEFT:
				pos = (pos-1+len(scenes))%len(scenes)
				#screen.fill((0,0,0)) 
			elif event.key == K_RIGHT:
				pos = (pos+1)%len(scenes)
				#screen.fill((0,0,0)) 
				
			