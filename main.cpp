#include <stdio.h>
#include <SDL/SDL.h>
#include <glm/glm.hpp>
#include <thread>
#include <mutex>
#include <vector>

#define WIDTH 512
#define HEIGHT 512
#define _STR(x) #x
#define STR(x) _STR(x)

using namespace glm;

// Pixel value
struct Pixel
{
	uint8_t r, g, b, a;
};

// A ray structure
struct Ray
{
	vec3 origin, direction;
};

// Material
struct Material
{
	// Surface color
	vec3 color;
	
	// Is the surface emissive
	bool emissive;
	
	// 0.0 => The surface is diffuse
	// 1.0 => The surface is a perfect mirror
	float shininess;
};

// Base class for geometry
class Geometry
{
public:
	virtual bool intersect(const Ray &ray, float &t0, vec3 &normal) const = 0;
};

// Sphere geometry
class Sphere : 
	public Geometry
{
public:
	vec3 center;
	float radius;
	float radius2;
	
	Sphere(const vec3 &center, float radius) : 
		center(center), radius(radius), radius2(radius*radius)
	{
		
	}
	
	virtual bool intersect(const Ray &ray, float &t0, vec3 &normal) const override
	{
		float t1; // solutions for t if the ray intersects

		// geometric solution
		vec3 L = this->center - ray.origin;
		float tca = dot(L, ray.direction);
		if (tca < 0) return false;
		float d2 = dot(L, L) - tca * tca;
		if (d2 > this->radius2) return false;
		float thc = sqrtf(this->radius2 - d2);
		t0 = tca - thc;
		t1 = tca + thc;
		
		normal = normalize(ray.origin + t0 * ray.direction - this->center);
		
		return true;
	}
};

// Plane geometry 
class Plane : 
	public Geometry
{
public:
	vec3 pos, normal;
	
	Plane(const vec3 &pos, const vec3 &normal) : 
		pos(pos), normal(normal)
	{
		
	}
	
	virtual bool intersect(const Ray &ray, float &t0, vec3 &normal) const override
	{
		float a = dot(ray.direction, this->normal);
		if(a == 0.0f)
			return false;
		t0 = dot(this->pos - ray.origin, this->normal) / a;
		normal = this->normal;
		return t0 >= 0;
	}
};

// Scene object
struct Object
{
	Material material;
	Geometry *geometry;
};

// The result screen
Pixel screen[WIDTH * HEIGHT];
// The scene
std::vector<Object> scene;

// Display the resulting image.
void display()
{
	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_Surface *surface = SDL_SetVideoMode(WIDTH, HEIGHT, 32, 0);
	
	memcpy(surface->pixels, screen, sizeof(screen));
	
	SDL_WM_SetCaption("Lylios Path Tracer", nullptr);
	
	SDL_Flip(surface);
	
	SDL_Event e;
	while(true)
	{
		if(SDL_PollEvent(&e))
		{
			if(e.type == SDL_QUIT)
				break;
			if((e.type == SDL_KEYDOWN) && (e.key.keysym.sym == SDLK_ESCAPE))
				break;
		}
		SDL_Delay(10);
	}
	
	SDL_Quit();
}

// Cast a ray into the scene
bool doRay(const Ray &ray, vec3 &hit, vec3 &normal, Material &surface)
{
	float minDist = FLT_MAX;
	
	size_t len = scene.size();
	for(size_t i = 0; i < len; i++)
	{
		float distance; vec3 n;
		const Object &object = scene[i];
		
		if(object.geometry->intersect(ray, distance, n) && (distance < minDist))
		{
			minDist = distance;
			surface = object.material;
			hit = ray.origin + ray.direction * distance;
			normal = n;
		}
	}
	
	if(minDist == FLT_MAX)
		return false;
		
	return true;
}

// Simple random generator
class Random
{
private:
	size_t m_z;
	size_t m_w;
public:
	Random()
	{
		this->seed(0);
	}
	
	void seed(uint32_t seed)
	{
		this->m_z = 0xAD9B6A3C ^ (seed);
		this->m_w = 0xCDCDBBDD ^ (seed >> 16);
	}
	
	// Generate uniform random in (0;1)
	double randf()
	{
		m_z = 36969 * (m_z & 65535) + (m_z >> 16);
		m_w = 18000 * (m_w & 65535) + (m_w >> 16);
		uint32_t u = (m_z << 16) + m_w;
		
		// The magic number below is 1/(2^32 + 2).
		// The result is strictly between 0 and 1.
		return (u + 1.0) * 2.328306435454494e-10;
	}	
};

#define USE_COSINE_SAMPLING

inline vec3 hemisphere(Random &rng, const vec3 &normal)
{
	float u1 = rng.randf();
	float u2 = rng.randf();
	
	// Use cosine sampling instead of linear sampling so we don't need to use the cosine weighting in the render equation
#if defined(USE_COSINE_SAMPLING)
	float r = sqrtf(u1);
	float theta = 2.0f * M_PI * u2;
	float x = r * cosf(theta);
	float y = r * sinf(theta);
	
	vec3 dir(x, y, sqrt(1 - u1));
#else
	// Use simple linear sampling
	float r = sqrtf(1.0f - u1 * u1);
    float phi = 2 * M_PI * u2;
    vec3 dir(cosf(phi) * r, sinf(phi) * r, u1);
#endif
///*
	// Return vector in the hemisphere of the normal 
	vec3 tangent;
	if(dot(normal, vec3(1,0,0)) > 0.8)
		tangent = cross(normal, vec3(1, 0, 0));
	else
		tangent = cross(normal, vec3(0, 0, 1));
	vec3 bitangent = cross(tangent, normal);
	
	return normalize(dir.x * tangent + dir.y * bitangent + dir.z * normal);

//*/
/*
    if(dot(dir, normal) < 0.0f)
		return -dir;
	else
		return dir;
//*/
}

// Focal length of the camera
const float focalLength = 1.2f;

// Screen aspect
const float aspect = (float)WIDTH / (float)HEIGHT;

// Number of samples for each pixel
const size_t numSamples = 500;

// Number of maximum rebounces
const size_t depth = 20;

// Camera setup
vec3 cPosition(2, 3, -5);
vec3 cForward(normalize(-cPosition)); // Just focus the camera on the origin
vec3 cRight = normalize(cross(cForward, vec3(0,1,0))); // Construct screen x axis 
vec3 cUp = normalize(cross(cForward, cRight)); 		   // Construct screen y axis

void render_scanline(Random &rng, size_t scanline)
{
	for(size_t x = 0; x < WIDTH; x++)
	{
		float fx = aspect * 2.0f * ((float)x / (float)WIDTH - 0.5f);
		float fy = 2.0f * (float)scanline / (float)(HEIGHT-1) - 1.0f;
		
		// Use 64 bit pixel values
		size_t color[] = { 0, 0, 0 };
		
		for(size_t sampleID = 0; sampleID < numSamples; sampleID++)
		{	
			// Construct ray from camera position
			Ray ray = {
				cPosition,
				normalize(focalLength * cForward + fx * cRight + fy * cUp)
			};
			vec3 sample(1.0f, 1.0f, 1.0f);
			for(size_t i = 0; i < depth; i++)
			{			
				vec3 hit, normal;
				Material surface;
				
				// Actual ray casting is here
				if(doRay(ray, hit, normal, surface) == false)
				{
					// Downscale the sample bei 0.3 for ambient lighting
					sample *= 0.3f;
					
					// If we haven't hit anything yet, set the result to black
					if(i == 0) sample *= 0.0f;
					break;
				}
				
				// Insert the hit surface into our resultig color
				sample *= surface.color;
				
				// Lights don't reflect light but emit, so quit here
				if(surface.emissive)
					break;

				// Path tracer core (reflect the ray in a direction between perfect reflection and perfect diffusion)
				ray.direction = mix(
					hemisphere(rng, normal), // Diffusion
					reflect(ray.direction, normal), // Reflection
					surface.shininess);
					
				// Place the ray to the hit position and move it further a bit so we don't collide with the hitpoint itself
				ray.origin = hit + 0.001f * ray.direction;
			}
			
			// Accumulate resulting sample values
			color[0] += 255 * sample.r;
			color[1] += 255 * sample.g;
			color[2] += 255 * sample.b;
		}
		
		// Write the pixel
		screen[WIDTH * scanline + x] = {
			(uint8_t)clamp<size_t>(color[2] / numSamples, 0, 255),
			(uint8_t)clamp<size_t>(color[1] / numSamples, 0, 255),
			(uint8_t)clamp<size_t>(color[0] / numSamples, 0, 255),
			255
		};
	}
}

// Mutex for synchronized line increasing
std::mutex mtx;

void render_scanlines(volatile size_t *current)
{
	bool quitafter = false;
	Random rng;
	do
	{
		// Increase the current scanline thread safe
		size_t scanline;
		{
			std::lock_guard<std::mutex> lock(mtx);
			scanline = *current;
			(*current)++;
			if((*current) >= HEIGHT)
				quitafter = true; // Check if we reached the end
		}
		if(quitafter == false)
		{
			//printf("Render Scanline %d...\n", (int)scanline);
			render_scanline(rng, scanline);
		}
	} while(quitafter == false);
	// Debug message for the random generator
	if(rng.randf() == rng.randf()) printf("rng broken! \n");
}

void render()
{
	// Start 4 rendering threads.
	
	volatile size_t currentScanline = 0;
	std::thread t1(render_scanlines, &currentScanline);
	std::thread t2(render_scanlines, &currentScanline);
	std::thread t3(render_scanlines, &currentScanline);
	std::thread t4(render_scanlines, &currentScanline);

	t1.join();
	t2.join();
	t3.join();
	t4.join();
}

// Built the scene
void genScene()
{
	// Ground
	scene.push_back(Object {
		Material {
			vec3(0.8f, 0.8f, 0.8f), false, 0.0f
		},
		new Plane(vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f))
	});
	
	// Lighted central sphere
	scene.push_back(Object {
		Material {
			10.0f * vec3(1.0f, 1.0f, 1.0f), true, 0.0f
		},
		new Sphere(vec3(0.0f, 1.5f, 0.0f), 0.4f)
	});
	
	// Red diffuse sphere
	scene.push_back(Object {
		Material {
			vec3(1.0f, 0.0f, 0.0f), false, 0.0f
		},
		new Sphere(vec3(2.5f, 0.8f, 0.0f), 1.0f)
	});
	
	// Green reflecting sphere
	scene.push_back(Object {
		Material {
			vec3(0.0f, 1.0f, 0.0f), false, 1.0f
		},
		new Sphere(vec3(-2.5f, 0.8f, 0.0f), 1.0f)
	});
}

int main(int argc, char **argv)
{
	genScene();
	
	render();

	display();
	
	return 0;
}

