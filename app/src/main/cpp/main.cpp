#include <iostream>
#include <array>
#include "raylib.h"
#include <algorithm>
#include <tuple>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"


class FluidPlane {
public:
    int size;
    float dt;
    float diff;
    float visc;

    float* temp_array;
    float* density;

    float* Vx;
    float* Vy;

    float* Vx0; 
    float* Vy0;

    FluidPlane(int size, float diff, float visc, float dt) {
        this->size = size;
        this->diff = diff;
        this->visc = visc;
        this->dt = dt;

        this->temp_array = new float[(size + 2) * (size + 2)] {};
        this->density = new float[(size + 2) * (size + 2)] {};

        this->Vx = new float[(size + 2) * (size + 2)] {};
        this->Vy = new float[(size + 2) * (size + 2)] {};

        this->Vx0 = new float[(size + 2) * (size + 2)] {};
        this->Vy0 = new float[(size + 2) * (size + 2)] {};

    }

    ~FluidPlane() {
        delete[] temp_array;
        delete[] density;
        delete[] Vx;
        delete[] Vy;
        delete[] Vx0;
        delete[] Vy0;
    }

    inline int IX(int i, int j) const {
        return i + (size + 2) * j;
    }

    void AddDensity(int x, int y, float amount) {
        int index = IX(x, y);
        this->density[index] += amount;

        if (this->density[index] >= 3000) this->density[index] = 3000;
        if (this->density[index] <= 0) this->density[index] = 0;
    }

    void AddVelocity(int x, int y, float amountX, float amountY) {
        int index = IX(x, y);
        this->Vx[index] += amountX;
        this->Vy[index] += amountY;
    }

    //void LinSolve(..) { .. }

    void SimulationStep(float dt) {
        int N = this->size;
        float visc = this->visc;
        float diff = this->diff;
        //float dt = this->dt;
        float* Vx = this->Vx;
        float* Vy = this->Vy;
     
        float* Vx0 = this->Vx0;
        float* Vy0 = this->Vy0;
  
        float* temp_array = this->temp_array;
        float* density = this->density;

        Diffuse(1, Vx0, Vx, visc, dt, N);
        Diffuse(2, Vy0, Vy, visc, dt, N);

        Project(Vx0, Vy0, Vx, Vy, N);

        Advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
        Advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

        Project(Vx, Vy, Vx0, Vy0, N);

        Diffuse(0, temp_array, density, diff, dt, N);
        Advect(0, density, temp_array, Vx, Vy, dt, N);

    }

    void SetBounds(int b, float* x, int N) {

        for (int i = 1; i <= N; i++) {
            x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
        }

        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    void Diffuse(int b, float* x, float* x0, float diff, float dt, int N){
        
        float a = dt * diff * N * N;

        //Gauss-Seidel Relaxation

        for (int k = 0; k < 20; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                        x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
                }
            }
            SetBounds(b, x, N);
        }

    }

    void Project(float* Vx, float* Vy, float* p, float* div, int N){
        int i, j, k;
        float h = 1.0 / N;

        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                div[IX(i, j)] = -0.5 * h * (Vx[IX(i + 1, j)] - Vx[IX(i - 1, j)] +
                    Vy[IX(i, j + 1)] - Vy[IX(i, j - 1)]);
                p[IX(i, j)] = 0;
            }
        }

        SetBounds(0, div, N); 
        SetBounds(0, p, N);

        // this might go to the function LinSolve(..){..}
        for (k = 0; k < 20; k++) {
            for (i = 1; i <= N; i++) {
                for (j = 1; j <= N; j++) {
                    p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                        p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
                }
            }
            SetBounds(0, p, N);
        }

        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                Vx[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
                Vy[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
            }
        }
        SetBounds(1, Vx, N);
        SetBounds(2, Vy, N);
    }
    
    void Advect(int b, float* d, float* d0, float* u, float* v, float dt, int N){

        int i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1;
        float dt0 = dt * N;

        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {

                x = i - dt0 * u[IX(i, j)]; 
                y = j - dt0 * v[IX(i, j)];

                if (x < 0.5) x = 0.5; 
                if (x > N + 0.5) x = N + 0.5; 
                
                i0 = (int)x; i1 = i0 + 1;
                
                if (y < 0.5) y = 0.5; 
                if (y > N + 0.5) y = N + 0.5; 
                
                j0 = (int)y; j1 = j0 + 1;

                s1 = x - i0; s0 = 1 - s1; 
                t1 = y - j0; t0 = 1 - t1;

                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                    s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }

        SetBounds(b, d, N);
    
    }


};


Color LerpColor(Color c1, Color c2, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    return Color{
        static_cast<unsigned char>((1 - t) * c1.r + t * c2.r),
        static_cast<unsigned char>((1 - t) * c1.g + t * c2.g),
        static_cast<unsigned char>((1 - t) * c1.b + t * c2.b),
        255
    };
}

void DrawCell(int i, int j, float val, int cellSize, int offset) {
    float maxDensity = 1;
    val = std::clamp(val / maxDensity, 0.0f, 1.0f);

    val *= 10;

    Color colorStart = { 212, 234, 98, 255 };
    Color colorMid = { 237, 113, 30, 255 };
    Color colorEnd = { 233, 41, 92, 255 };

    Color col;
    if (val < 0.5f) {
        float t = val * 2.0f;
        col = LerpColor(colorStart, colorMid, t);
    }
    else {
        float t = (val - 0.5f) * 2.0f;
        col = LerpColor(colorMid, colorEnd, t);
    }

    DrawRectangle(i * cellSize + offset, j * cellSize, cellSize, cellSize, col);
}


void DrawSimulationCells(FluidPlane& fluidPlane, int width, int height, float fading, int offset) {

    int size = fluidPlane.size;
    float rectSide = width / float(size);


    for (int i = 1; i <= size; i++) {
        for (int j = 1; j <= size; j++) {

            const auto val = fluidPlane.density[fluidPlane.IX(i, j)];
            //DrawRectangleLines((i - 1) * rectSide, (j - 1) * rectSide, rectSide, rectSide, RED);
            DrawCell((i - 1), (j - 1), val, rectSide, offset);

            fluidPlane.AddDensity(i, j, -fading);
            
        }
    }

}

std::tuple<int, int> GetCellFromPosition(float x, float y, int size, int screenWidth, int screenHeight, int offset) {
    float rectSide = screenWidth / float(size);

    int i = (x - offset) / rectSide;
    int j = y / rectSide;

    return std::tuple<int, int>{i, j};
}

int main(){

    InitWindow(0, 0, "FluidSim");

    int screenWidth = GetScreenWidth();
    int screenHeight = GetScreenHeight();
    
    int fluidPlaneWidth = screenHeight;
    int fluidPlaneHeight = screenHeight;
    int offset = screenWidth - fluidPlaneWidth;

    int FPS = 60;
    SetTargetFPS(FPS);

    Rectangle fluidBounds = Rectangle(0 + offset, 0, fluidPlaneHeight, fluidPlaneWidth);

    int size = 100;
    
    float diff = 0.0001;
    float diffLog = -4;

    float viscosity = 0.0001;
    float viscosityLog = -4;
    
    float fadingRate = 0.000001;
    float fadingLog = -6;
    
    float density = 10.0;

    float dt = 1.0 / FPS;

    FluidPlane fluid = FluidPlane(size, diff, viscosity, dt);

    //GuiSetStyle(DEFAULT, TEXT_SIZE, 12);
    
    int gui_x = offset * 0.3; //120;
    int gui_w = offset * 0.5; //200;
    int gui_h = screenHeight * 0.05; //40 
    int gui_spacing = screenHeight * 0.04;


    Vector2 oldMousePos = GetMousePosition();

    while (!WindowShouldClose()) {
        //dt = GetFrameTime();
        std::cout << GetFPS() << std::endl;

        BeginDrawing();

        ClearBackground(BLACK);
        
        fluid.SimulationStep(dt);
        DrawSimulationCells(fluid, fluidPlaneWidth, fluidPlaneHeight, fadingRate, offset);

        Vector2 newMousePos = GetMousePosition();

        // the distance our cursor has moved since the last frame
        float distanceX = newMousePos.x - oldMousePos.x;
        float distanceY = newMousePos.y - oldMousePos.y;

        const auto [i, j] = GetCellFromPosition(newMousePos.x, newMousePos.y, size, fluidPlaneWidth, fluidPlaneHeight, offset);

        if (CheckCollisionPointRec(newMousePos, fluidBounds)){
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
                fluid.AddDensity(i, j, density);
                fluid.AddVelocity(i, j, distanceX, distanceY);
            }
            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
                fluid.AddVelocity(i, j, distanceX, distanceY);
            }
        // Update mouse's previous position
        oldMousePos = newMousePos;
        }

        GuiSliderBar(Rectangle(gui_x, gui_h + gui_spacing, gui_w, gui_h), "Diffusion Rate", TextFormat("%.6f", diff), &diffLog, -6.0f, -1.0f);
        diff = powf(10.0f, diffLog);
        fluid.diff = diff;
        
        GuiSliderBar(Rectangle(gui_x, 2*(gui_h + gui_spacing), gui_w, gui_h), "Viscosity", TextFormat("%.6f", viscosity), &viscosityLog, -6.0f, -1.0f);
        viscosity = powf(10.0f, viscosityLog);
        fluid.visc = viscosity;

        GuiSliderBar(Rectangle(gui_x, 3 * (gui_h + gui_spacing), gui_w, gui_h), "Fluid Density", TextFormat("%.1f", density), &density, 1.0f, 10.0f);
        
        GuiSliderBar(Rectangle(gui_x, 4 * (gui_h + gui_spacing), gui_w, gui_h), "Fading Rate", TextFormat("%.6f", fadingRate), &fadingLog, -6.0f, -1.0f);
        fadingRate = powf(10.0f, fadingLog);

        EndDrawing();

    }

    CloseWindow();
    return 0;


}
