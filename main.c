#ifdef _WIN32 // грубый способ определения; на деле разница не в ОС, а в компиляторe
	#include <SDL.h>
	#include <SDL_image.h>
	#include <SDL_timer.h>
#else
	#include <SDL2/SDL.h>
	#include <SDL2/SDL_image.h>
	#include <SDL2/SDL_timer.h>
#endif

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <string.h>

/* TILESIZE (и SQRTTILE) лучше не трогать - или быть 
готовым гадать, где 32 заменить, а где не надо */
/* При изменении размера окна не забыть поменять 
скайбокс (или добавить масштабируемость) */

#define PI 3.1415926535 // значение числа пи
#define DEG 0.01745329	// знаечние 1 градуса в радианах

#include "./var.conf"
#include "levels/level_test.dat"

#define TILESIZE 32		// размер одной клетки карты; также "реальная" высота стены	
#define SQRTTILE 5		// двоичный логарифм TILESIZE - используется для побитовых операций

#if (MAPX > MAPY)		// дальность прорисовки
	#define DOF MAPX
#else
	#define DOF MAPY
#endif

// ZWALL - noclip
// ZINFRA - отображение всех спрайтов на карте
int cheat_active[2] = {0,0,};
char cheatcode[2][8] = {"wall", "infra"};
char cheat[8];
int cheatlen = 0;

time_t systime;

typedef struct { float x, y; } vec;
vec dir = {TILESIZE, 0};					// вектор направления
vec plane = {0, TILESIZE*tan(FOV/2*DEG)};	// вектор пл-ти камеры
// double invDet = 1.0 / (plane.x * dir.y - plane.y * dir.x); 	// обратный определитель
double invDet = 1.0 / (- TILESIZE*TILESIZE*tan(FOV/2*DEG));		// для отрисовки спрайтов

float dist(float x1, float y1, float x2, float y2)
{
	return ( sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) );
}

void clearCheatCache()
{
	memset(cheat, 0, sizeof(cheat));
	cheat[0] = '\0';
	cheatlen = 0;
	SDL_StopTextInput();
}

Uint32 mixSDLColors(Uint32 back, Uint32 front, SDL_PixelFormat *fmt, int alpha)
{
	if(alpha > 127)
		return front;

	Uint32 ch1, ch2, rmask, gmask, bmask, sum = 0;
	rmask = fmt->Rmask;
	gmask = fmt->Gmask;
	bmask = fmt->Bmask;

	ch1 = back & rmask;
	ch2 = front & rmask;
	sum |= (ch1*(128-alpha) + ch2*alpha) >> 7 & rmask;

	ch1 = back & gmask;
	ch2 = front & gmask;
	sum |= (ch1*(128-alpha) + ch2*alpha) >> 7 & gmask;

	ch1 = back & bmask;
	ch2 = front & bmask;
	sum |= (ch1*(128-alpha) + ch2*alpha) >> 7 & bmask;

	sum |= front & !(rmask | gmask | bmask);
	return sum;
}

void rotateVector(vec *v, float ang)
{
	/* Матрица поворота в левосторонней системе координат: */
	// |  cos(a)  sin(a) | 
	// | -sin(a)  cos(a) |
	float xr = (v->x)*cos(ang) - (v->y)*sin(ang);
	float yr = (v->x)*sin(ang) + (v->y)*cos(ang);
	v->x = xr;
	v->y = yr;
}

void sortSprites(int order[], float dist[], int amount)
{
	double factor = 1.24733095;
	int step = amount-1;
	double tmp;
	while(step >= 1)
	{
		for(int i = 0; i+step < amount; i++)
		{
			if(dist[i] < dist[i+step])
			{
				tmp = dist[i];
				dist[i] = dist[i+step];
				dist[i+step] = tmp;

				tmp = (double)order[i];
				order[i] = order[i+step];
				order[i+step] = (int)tmp;
			}
		}
		step /= factor;
	}
}

void scalePPM(int *img)
{
	;
}

int processEvents(SDL_Window *window, struct player *P, int map[], int *mpbool, clock_t *prevframe)
{
	// общеигровое время unix epoch
	systime = time(NULL);
	// Вычисление FPS
	clock_t frame;
	float fps;
	
	frame = clock();
	fps = CLOCKS_PER_SEC/(double)(frame - *prevframe);
	*prevframe = frame;
	printf("%d\t\t\r", (int)fps);

	int done = 0;
	SDL_Event event;

	int xoff = 0;
	if((P->dx)>0) {xoff = 20;} else {xoff = -20;}

	int yoff = 0;
	if((P->dy)>0) {yoff = 20;} else {yoff = -20;}

	int px64 = (int)(P->x) >> SQRTTILE;
	int py64 = (int)(P->y) >> SQRTTILE;

	int pxadd = (int)(P->x + xoff) >> SQRTTILE, pyadd = (int)(P->y + yoff) >> SQRTTILE;
	int pxsub = (int)(P->x - xoff) >> SQRTTILE, pysub = (int)(P->y - yoff) >> SQRTTILE;

	while(SDL_PollEvent(&event))
	{
		switch(event.type)
		{	
			case SDL_WINDOWEVENT_CLOSE:
			{
				if(window)
				{
					SDL_DestroyWindow(window);
					window = NULL;
					done = 1;
				}
			}
			break;
			
			case SDL_KEYDOWN:
			{
				switch(event.key.keysym.sym)
				{
					case SDLK_ESCAPE:
						done = 1;
						break;
					case SDLK_m: // карта
						*mpbool = !*mpbool;
						break;
					case SDLK_e: // взаимодействие
						int xeoff = 0; if((P->th)<(PI/4) || (P->th)>(7*PI/4))  {xeoff = 25;} 
								  else if((P->th)>(3*PI/4) && (P->th)<(5*PI/4)){xeoff = -25;}

						int yeoff = 0; if((P->th)>(PI/4) && (P->th)<(3*PI/4))   {yeoff = 25;} 
								  else if((P->th)>(5*PI/4) && (P->th)<(7*PI/4)) {yeoff = -25;}
					
						int exadd = ((int)P->x + xeoff) >> SQRTTILE, eyadd = ((int)P->y + yeoff) >> SQRTTILE;

						if(map[eyadd*MAPX + exadd] == 4) // дверь
						   map[eyadd*MAPX + exadd] = 0;
						break;
					case SDLK_z:
						clearCheatCache();
						SDL_StartTextInput();
				}
			}
			break;

			case SDL_TEXTINPUT:
				strcat(cheat, event.text.text);
				for(int i = 0; i<2; i++)
					if(strcmp(cheat, cheatcode[i])==0)
					{
						cheat_active[i] = !cheat_active[i];
						if(cheat_active[i])
							printf("z%s enabled\n", cheatcode[i]);
						else
							printf("z%s disabled\n", cheatcode[i]);
						clearCheatCache();
					}

				cheatlen++;
				if(cheatlen >= 8) 
				{
					clearCheatCache();
				}
				break;
			
			case SDL_QUIT:
				done = 1;
			break;
		}
	}

	const Uint8 *state = SDL_GetKeyboardState(NULL);

	if(state[SDL_SCANCODE_W]) {
		if(map[py64*MAPX + pxadd]==0 || cheat_active[0]) {P->x += P->dx;}
		if(map[pyadd*MAPX + px64]==0 || cheat_active[0]) {P->y += P->dy;}
	}

	if(state[SDL_SCANCODE_S]) {
		if(map[py64*MAPX + pxsub]==0 || cheat_active[0]) {P->x -= P->dx;}
		if(map[pysub*MAPX + px64]==0 || cheat_active[0]) {P->y -= P->dy;}
	} 

	if(state[SDL_SCANCODE_A]) {
		P->th -= (RADSPD/fps);
		rotateVector(&dir,   -RADSPD/fps);
		rotateVector(&plane, -RADSPD/fps);
	}

	if(state[SDL_SCANCODE_D]) {
		P->th += (RADSPD/fps);
		rotateVector(&dir,   RADSPD/fps);
		rotateVector(&plane, RADSPD/fps); 
	}

	if(P->th < 0)
		P->th = 2*PI - 0.0001;
	if(P->th > 2*PI)
		P->th = 0.0;

	P->dx = cos(P->th)*(SPEED/fps);
	P->dy = sin(P->th)*(SPEED/fps);

	return done;
}

void drawMinimap(SDL_Renderer *renderer, int map[], struct player P)
{
	for(int y = 0; y<MAPY; y++)
		for(int x = 0; x<MAPX; x++)
		{
			if(map[y*MAPX + x] > 0)
			{
				SDL_SetRenderDrawColor(renderer, 255, 255, 255, 128);
			} else {
				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 128);
			}

				SDL_Rect wall = {x*TILESIZE/MAPSCALE, y*TILESIZE/MAPSCALE, (TILESIZE-4)/MAPSCALE, (TILESIZE-4)/MAPSCALE};
				SDL_RenderFillRect(renderer, &wall);
		}

	SDL_Vertex arrow[3];
	arrow[0].position.x = (P.x+7*cos(P.th-PI/2))/MAPSCALE;
	arrow[0].position.y = (P.y+7*sin(P.th-PI/2))/MAPSCALE;
	arrow[0].color.r = 255;
	arrow[0].color.g = 220;
	arrow[0].color.b = 0;
	arrow[0].color.a = 128;

	arrow[1].position.x = (P.x+7*cos(P.th+PI/2))/MAPSCALE;
	arrow[1].position.y = (P.y+7*sin(P.th+PI/2))/MAPSCALE;
	arrow[1].color.r = 255;
	arrow[1].color.g = 220;
	arrow[1].color.b = 0;
	arrow[1].color.a = 128;

	arrow[2].position.x = (P.x+20*cos(P.th))/MAPSCALE;
	arrow[2].position.y = (P.y+20*sin(P.th))/MAPSCALE;
	arrow[2].color.r = 255;
	arrow[2].color.g = 220;
	arrow[2].color.b = 0;
	arrow[2].color.a = 128;

	if(cheat_active[1])
	{
		SDL_SetRenderDrawColor(renderer, 255, 0, 0, 128);
		for(int i = 0; i<sizeof(sp)/sizeof(sp[0]); i++)
		{
			SDL_Rect entity = {sp[i].x/MAPSCALE-4, sp[i].y/MAPSCALE-4, 4, 4};
			SDL_RenderFillRect(renderer, &entity);
		}
	}

	SDL_RenderGeometry(renderer, NULL, arrow, 3, NULL, 0);		
}

void drawRays(SDL_Renderer *renderer, void *tpixels, SDL_Surface *surface, struct player P, int map[], int mapF[], int mapC[], float zbuffer[])
{
	int mx, my, mp, dof, texture_iv = 4, texture_ih = 4, texture_i;
	float rx, ry, ra, xo, yo, px, py, hlen, vlen, sinra, cosra, 
		  hx = 0, hy = 0, vx = 0, vy = 0, shade_h = 1, shade_v = 1, rdist;

	Uint32 *pixels = (Uint32 *)malloc(WIDTH*HEIGHT*sizeof(Uint32));

	ra=P.th-FOV/2*DEG;
	if(ra<0)
		ra += 2*PI;
	if(ra>2*PI)
		ra -= 2*PI;

	px = P.x;
	py = P.y;
	float antan;
	for(int r = 0; r<WIDTH/PIXELSIZE; r++)
	{
		// ----- по горизонтальным -----
		hlen = 100000;
		dof = 0;
		sinra = sin(ra);
		antan=-1.0/tan(ra);
		if(sinra<-0.001)	// "вверх"
		{
			ry = (((int)py>>SQRTTILE)<<SQRTTILE)-0.0001;
			rx = (py-ry)*antan+px;
			yo = -TILESIZE;
			xo = -yo*antan;
		} else {
		if(sinra>0.001)
		{
			ry = (((int)py>>SQRTTILE)<<SQRTTILE)+TILESIZE;
			rx = (py-ry)*antan+px;
			yo = TILESIZE;
			xo = -yo*antan;
		}
		else
		{
			rx = px;
			ry = py;
			dof = DOF;
		}}

		while(dof<DOF)
		{	
			mx = (int)(rx)>>SQRTTILE;
			my = (int)(ry)>>SQRTTILE;
			mp = my*MAPX+mx;

			if(mp > 0 && mp < MAPX*MAPY && map[mp]>0){
				hx = rx;
				hy = ry;
				hlen = dist(px, py, hx, hy);
				texture_ih = map[mp]-1;
				shade_h = shades[mp];
				break;
			} else {
				rx += xo;
				ry += yo;
				dof++;
			}
		}
		// ----- по вертикальным -----
		dof = 0;
		vlen = 100000;
		cosra = cos(ra);
		antan=-tan(ra);
		if(cosra<-0.001)	// "влево"
		{
			rx = (((int)px>>SQRTTILE)<<SQRTTILE)-0.0001;
			ry = (px-rx)*antan+py;
			xo = -TILESIZE;
			yo = -xo*antan;
		} else {
		if(cosra>0.001)
		{
			rx = (((int)px>>SQRTTILE)<<SQRTTILE)+TILESIZE;
			ry = (px-rx)*antan+py;
			xo = TILESIZE;
			yo = -xo*antan;
		}
		else
		{
			rx = px;
			ry = py;
			dof = DOF;
		}}

		while(dof<DOF)
		{	
			mx = (int)(rx)>>SQRTTILE;
			my = (int)(ry)>>SQRTTILE;
			mp = my*MAPX+mx;

			if(mp > 0 && mp < MAPX*MAPY && map[mp]>0){
				vx = rx;
				vy = ry;
				vlen = dist(px, py, vx, vy);
				texture_iv = map[mp]-1;
				shade_v = shades[mp];
				break;
			} else {
				rx += xo;
				ry += yo;
				dof++;
			}
		}

		float shade = 1;
		if(vlen < hlen) 
		{ 
			rx = vx; ry = vy; rdist = vlen; texture_i = texture_iv; shade = shade_v;
		}
		else { 
			rx = hx; ry = hy; rdist = hlen; texture_i = texture_ih; shade = shade_h*0.6;
		}

		// ----- Рисуем стены -----

		float corran = P.th - ra;
		if(corran<0)
			corran += 2*PI;
		if(corran>2*PI)
			corran -= 2*PI;
		rdist *= cos(corran); 						// рыбий глаз

		float lineh = (TILESIZE*HEIGHT)/rdist;
		float ty_step = 32.0 / lineh;
		float ty_offset = 0;
		
		if(lineh>HEIGHT){
			ty_offset = (lineh-HEIGHT)/2.0;
			lineh=HEIGHT;
		}
		float lineo = HEIGHT/2 - lineh/2;
		
		float ty = ty_offset*ty_step; 				// Текстурирование
		float tx = 0;
		if(vlen > hlen) {
			tx = (int)(rx/(TILESIZE/32)) % 32;
			if(ra < 180*DEG) {tx = 31-tx;}
		} else {
			tx = (int)(ry/(TILESIZE/32)) % 32;
			if(ra>90*DEG && ra < 270*DEG) {tx = 31-tx;}
		}

		int y = 0;
		for(y = lineo; y<lineo+lineh; y++, ty+=ty_step)
		{
			int index = ( ((int)ty*32 + (int)tx)*3+(texture_i*32*32*3) )%(32*192*3);
			char red = textures[index]*shade;
			char green = textures[index+1]*shade;
			char blue = textures[index+2]*shade;
	
			for(int i = 0; i<PIXELSIZE; i++)
				pixels[y*WIDTH + r*PIXELSIZE +i] = SDL_MapRGB(surface->format, red, green, blue);
				/*pixels[y*WIDTH + r*PIXELSIZE +i] = mixSDLColors(SDL_MapRGB(surface->format, red, green, blue),
																0,
																surface->format,
																128*(1-lineh/HEIGHT));*/
		}
 
		//   ----- рисуем пол -----
		for(y = lineo+lineh; y<HEIGHT; y++)
		{
			float dy = y - (HEIGHT/2.0);
			tx = P.x  + cos(ra)*(0.4941406*HEIGHT)*32/(dy*cos(corran));
			ty = P.y  + sin(ra)*(0.4941406*HEIGHT)*32/(dy*cos(corran));
			int array_i = ( (int)(ty/(float)TILESIZE)*MAPX + (int)(tx/(float)TILESIZE) )%(MAPX*MAPY);
			int mp = mapF[array_i]*32*32;

			int index = (((int)(ty)&31)*32 + ((int)(tx)&31))*3 + mp*3;
			index %= sizeof(textures)/sizeof(textures[0]);							////
			char red = textures[index]*0.7*shades[array_i];
			char green = textures[index+1]*0.7*shades[array_i];
			char blue = textures[index+2]*0.7*shades[array_i];

			for(int i = 0; i<PIXELSIZE; i++)
				pixels[y*WIDTH + r*PIXELSIZE +i] = SDL_MapRGB(surface->format, red, green, blue);
				/*pixels[y*WIDTH + r*PIXELSIZE +i] = mixSDLColors(SDL_MapRGB(surface->format, red, green, blue),
																0,
																surface->format,
																128*(1-dy/(HEIGHT/2)));*/

			// ----- рисуем потолок -----
			mp = mapC[array_i];

			if(mp > 0)
			{
				index = (((int)(ty)&31)*32 + ((int)(tx)&31))*3 + (mp-1)*32*32*3;
				index %= sizeof(textures)/sizeof(textures[0]);							////
				red = textures[index]*0.7*shades[array_i];
				green = textures[index+1]*0.7*shades[array_i];
				blue = textures[index+2]*0.7*shades[array_i];

				for(int i = 0; i<PIXELSIZE; i++)
					pixels[(HEIGHT-y)*WIDTH + r*PIXELSIZE +i] = SDL_MapRGB(surface->format, red, green, blue);
					/*pixels[(HEIGHT-y)*WIDTH + r*PIXELSIZE +i] = mixSDLColors(SDL_MapRGB(surface->format, red, green, blue),
																0,
																surface->format,
																128*(1-dy/(HEIGHT/2)));*/
			} else {
				for(int x = r*PIXELSIZE; x < r*PIXELSIZE+PIXELSIZE; x++)
				{
					int xo = (int)(-P.th*WIDTH/(FOV*DEG)) - x;
					if(xo < 0) {xo += WIDTH;}
					xo %= WIDTH; 
					index = (y*WIDTH + xo)*3;
					red = sky[index];
					green = sky[index+1];
					blue = sky[index+2];
					pixels[(HEIGHT-y)*WIDTH + x] = SDL_MapRGB(surface->format, red, green, blue);
				}
			}
		}

		for(int p = 0; p < PIXELSIZE; p++){
			zbuffer[(int)(r*PIXELSIZE)+p] = rdist;
		}

		ra += (double)FOV/((double)WIDTH/(double)PIXELSIZE)*DEG;
		if(ra<0)
			ra += 2*PI;
		if(ra>2*PI)
			ra -= 2*PI;
	}

	memcpy(tpixels, pixels, (size_t)(surface->pitch * surface->h));
	free(pixels);
}

void drawSprites(SDL_Renderer *renderer, void *tpixels, SDL_Surface *surface, struct player P, float zbuffer[])
{
	Uint32 *pixels = (Uint32 *)malloc(WIDTH*HEIGHT*sizeof(Uint32));

	int spnum = sizeof(sp)/sizeof(sp[0]);
	int  sp_order[spnum];
	float sp_dist[spnum];
	for(int i = 0; i < spnum; i++)
	{
		sp_order[i] = i;
		sp_dist[i] = ( (P.x - sp[i].x)*(P.x - sp[i].x) + (P.y - sp[i].y)*(P.y - sp[i].y) );
	}
	sortSprites(sp_order, sp_dist, spnum);

	for(int i = 0; i < spnum; i++)
	{
		if(sp[sp_order[i]].state == 0)
			continue;
		double relX = sp[sp_order[i]].x - P.x;		// относительные координаты спрайта
		double relY = sp[sp_order[i]].y - P.y;
		int shadeindex = ( (int)(sp[sp_order[i]].y/(float)TILESIZE)*MAPX + (int)(sp[sp_order[i]].x/(float)TILESIZE) )%(MAPX*MAPY);
		/* трансформируем координаты в экранные засчет 
		   умножения на обратную матрицу камеры */
		//								| нпY	-нпX |
		// -1/(плX * нпY - нпX * плY) * |			 |
		//								| -плY	 плX |
		// где (напX, напY) - вектор направления, (плX, плY) - вектор плоскости камеры
		//  __     ___
		// |пл| = |нап| * tg(FOV/2)
		// |нап| = HEIGHT
		// нап = (HEIGHT, 0)
		// пл = (0, HEIGHT*tg(FOV/2))

		// double corr = 1.045;
		//								| нпY  -нпХ	|	| x_отн	|	| x_тр |
		// 1/(плХ * нпY - нпХ * плY) *  |			| * |		| = |	   |
		//								| плY   плX |	| y_отн |	| y_тр |
		// x_тр - координата X в плоскости камеры, y_тр - фактор масштабирования
		double trX = invDet * (dir.y * relX - dir.x * relY);	// трансформируем координаты в координаты камеры
		double trY = invDet * ((-plane.y) * relX + plane.x * relY);
		int vmove = (int)( sp[sp_order[i]].offset / trY);

		int screenX = (int)((WIDTH/2) * (1 + trX / trY));

		int sp_height = abs((int)(HEIGHT / trY)) / sp[sp_order[i]].scaleH;
		int startY = HEIGHT/2 - sp_height/2 + vmove; if(startY < 0)  {startY = 0;}
		int endY = HEIGHT/2 + sp_height/2 + vmove; if(endY >= HEIGHT){endY = HEIGHT-1;}

		int sp_width = abs((int)(WIDTH / trY)) / sp[sp_order[i]].scaleV; // возможно, тут WIDTH?
		int startX = screenX - sp_width/2; if(startX < 0) {startX = 0;}
		int endX = screenX + sp_width/2; if(endX >= WIDTH){endX = WIDTH-1;}

		for(int s = startX; s < endX; s++)
		{
			int tx = (int)((s - (screenX - sp_width/2)) * 64 / sp_width);
			
			if(trY > 0 && s > 0 && s < WIDTH && trY*TILESIZE < zbuffer[s])
			{
				for(int y = startY; y < endY; y++)
				{
					int d = ((y-vmove)*256 - HEIGHT*128 + sp_height*128);
					int ty = (d * 64 / sp_height)/256;
					
					int index = (((int)(ty)&63)*64 + ((int)(tx)&63))*3 + (systime%4)*64*64*3;
					int red = sprite[index];
					int green = sprite[index+1];
					int blue = sprite[index+2];

						if(!(red == 255 && green == 0 && blue == 255)) // невидимый цвет (255, 0, 255)
							/* pixels[y*WIDTH + s] = SDL_MapRGB(surface->format, red*shades[shadeindex],
																green*shades[shadeindex],
																blue*shades[shadeindex]);            */
							pixels[y*WIDTH + s] = mixSDLColors(pixels[y*WIDTH + s], 
															   SDL_MapRGB(surface->format, 
						 				   						 		  red *shades[shadeindex], 
										   			  					  green *shades[shadeindex], 
										   			  					  blue *shades[shadeindex]),
															   surface->format, 
															   sp[sp_order[i]].opacity);

				} 
			}
		}
	}
	memcpy(tpixels, pixels, (size_t)(surface->pitch * surface->h));
	free(pixels);
}

int main(int argc, char *argv[])
{
	SDL_Window *window = NULL;
	SDL_Renderer *renderer = NULL;
	SDL_Surface *surface = NULL;
	SDL_Texture *texture = NULL;
	
	SDL_Init(SDL_INIT_VIDEO);

	window = SDL_CreateWindow("Raycast",
				  SDL_WINDOWPOS_CENTERED,
				  SDL_WINDOWPOS_CENTERED,
				  WIDTH,
				  HEIGHT,
				  0);

	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

	surface = SDL_GetWindowSurface(window);
	// SDL_SetSurfaceRLE(surface, 1);
	texture = SDL_CreateTexture(renderer, SDL_GetWindowPixelFormat(window), SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);

	clock_t prevframe = clock();
	clock_t *pnt_pframe = &prevframe;

	int done = 0;	// bool завершения программы
	int mpbool = 0; // bool вкл/выкл карты
	
	void *tpixels = NULL;
	int pitch = 0;
	float zbuffer[WIDTH*PIXELSIZE+1] = {0};

	while(!done)
	{
	// --- ПРОВЕРКА СОБЫТИЙ ---
		if(processEvents(window, &P, layoutW, &mpbool, pnt_pframe) == 1)
			done = 1;
	
	// --- РЕНДЕР --- 
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // фон
		SDL_RenderClear(renderer);

		SDL_LockTexture(texture, NULL, &tpixels, &pitch);
		// drawSky(tpixels, surface, P);
		drawRays(renderer, tpixels, surface, P, layoutW, layoutF, layoutC, zbuffer);
		drawSprites(renderer, tpixels, surface, P, zbuffer);
		SDL_UnlockTexture(texture);

		SDL_RenderCopy(renderer, texture, NULL, NULL);
			
		if(mpbool)
			drawMinimap(renderer, layoutW, P);
	
		SDL_RenderPresent(renderer);
		// SDL_Delay(DELAY);
	}

	SDL_DestroyWindow(window);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyTexture(texture);

	SDL_Quit();
	printf("\n\n!%d!", errno);
	return 0;
}