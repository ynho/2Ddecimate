#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#define SIDE_1 1
#define SIDE_2 2
#define SIDE_3 4
#define SIDE_4 8
#define SIDE_5 16
#define SIDE_6 32

#define MX SIDE_1
#define PX SIDE_2
#define MY SIDE_3
#define PY SIDE_4
#define MZ SIDE_5
#define PZ SIDE_6


/*
for each cube c:
  for each adjacent face f of c:
    m <- build index mapping of f to 2D mesh
for each joint j (8 total):
  v <- extract 2D geometry
  s <- simplify v
for each cube c:
  for each adjacent face f of c:
    reconnect local vertices to simplified s
    reconnect adjacent cubes... ? kinda already done if you just merge concatenate all the arrays except shared vertices & indices
  end;
end;
*/

struct tri {
    float area;
    int a, b;
};

static int ran (int a, int b) {
    float r = (float)rand () / RAND_MAX;
    r *= (b - a);
    r += a;
    return r + 0.5;
}

static float heron (int x[2], int y[2], int z[2]) {
    float a, b, c;
    float q, r;
    q = x[0] - y[0], r = x[1] - y[1];
    a = sqrt (q * q + r * r);
    q = y[0] - z[0], r = y[1] - z[1];
    b = sqrt (q * q + r * r);
    q = z[0] - x[0], r = z[1] - x[1];
    c = sqrt (q * q + r * r);
    float s = (a + b + c) * 0.5;
    return sqrt (s * (s - a) * (s - b) * (s - c));
}

static inline void v_connect (struct tri *tris, int a, int b) {
    if (tris[a].a == -1)
        tris[a].a = b;
    else /* if (tris[a].b == -1) */
        tris[a].b = b;
    /* else non-manifold input */
}

static void compute_area (struct tri *tri, int self, int *vertices, int side) {
    int x, y;
    if (side & (MX | PX)) {
        x = 1, y = 2;
    } else if (side & (MY | PY)) {
        x = 0, y = 2;
    } else {
        x = 0, y = 1;
    }
    if (tri->a < 0 || tri->b < 0)
        tri->area = -1.0;
    else {
        int a[2], b[2], c[2];
        a[0] = vertices[self * 3 + x];
        a[1] = vertices[self * 3 + y];
        b[0] = vertices[tri->a * 3 + x];
        b[1] = vertices[tri->a * 3 + y];
        c[0] = vertices[tri->b * 3 + x];
        c[1] = vertices[tri->b * 3 + y];
        tri->area = heron (a, b, c);
    }
}

static void compute_areas (struct tri *origin, struct tri **tris, int n_tris,
                           int *vertices, int side) {
    for (int i = 0; i < n_tris; i++) {
        int index = (int)(tris[i] - origin);
        compute_area (tris[i], index, vertices, side);
    }
}

static int secondary_condition(int *vertices, int a, int b) {
    if (vertices[a * 3] == vertices[b * 3]) {
        if (vertices[a * 3 + 1] == vertices[b * 3 + 1])
            return vertices[a * 3 + 2] < vertices[b * 3 + 2];
        else
            return vertices[a * 3 + 1] < vertices[b * 3 + 1];
    } else
        return vertices[a * 3] < vertices[b * 3];
}

static void sort_areas (struct tri **tris, int *vertices, int n_tris) {
    struct tri *x = NULL;
    for (int i = 1; i < n_tris; i++) {
        int j = i;
        while ((tris[j]->area < tris[j - 1]->area ||
                (tris[j]->area == tris[j - 1]->area &&
                 secondary_condition(vertices, j, j - 1)))
               && j > 1) {
            x = tris[j];
            tris[j] = tris[j - 1];
            tris[j - 1] = x;
            j--;
        }
    }
}

static void shift_left (struct tri **tris, int from, int to) {
    for (int i = from; i < to; i++)
        tris[i] = tris[i + 1];
}

void reconnect (struct tri *t, int old, int new) {
    if (t->a == old)
        t->a = new;
    else // if (t->b == old)
        t->b = new;
}

void pa (struct tri **tris, int n) {
    printf ("areas : ");
    for (int i = 0; i < n; i++)
        printf ("%.2f\t", tris[i]->area);
    printf ("\n");
}

void decimate (int n_vertices, int n_indices, int *vertices, int *indices,
               int *merge, int side) {
    int i;
    /* so here, n_vertices is actually correct but since we are getting
       all the vertices from the whole cube, hacking in the (tweaked) indices
       number is better (smaller number) */
    /* struct tri *tris = malloc (n_vertices * sizeof *tris); */
    struct tri *tris = malloc (n_vertices * sizeof *tris);
    size_t n_sorted = n_indices / 2; /* devided by 2 because 1 edge has 2 vertices */
    struct tri **sorted = malloc (n_sorted * sizeof *tris);

    for (i = 0; i < n_vertices; i++) {
        tris[i].a = tris[i].b = -1;
    }
    for (i = 0; i < n_sorted; i++) {
        int a = indices[i * 2], b = indices[i * 2 + 1];
        v_connect (tris, a, b);
        v_connect (tris, b, a);
        sorted[i] = &tris[a];
    }

    /* compute all areas */
    compute_areas (tris, sorted, n_sorted, vertices, side);

    /* sort all areas */
    sort_areas (sorted, vertices, n_sorted);

    int reduces = n_indices / 4;

    for (i = 0; i < reduces; i++) {
        /* skip borders */
        int offset = 0;
        /* TODO: can overflow */
        while (sorted[offset]->area < 0.0)
            offset++;

        int origin = (int)(sorted[offset] - tris);

        /* reduce least significant vertex */
        int a = tris[origin].a, b = tris[origin].b;

        merge[origin] = a;
        /* tris[origin].area = -1.0; */

        reconnect (&tris[a], origin, b);
        reconnect (&tris[b], origin, a);

        compute_area (&tris[a], a, vertices, side);
        compute_area (&tris[b], b, vertices, side);

        n_sorted--;
        shift_left (sorted, offset, n_sorted);
        sort_areas (&sorted[offset], vertices, n_sorted - offset);
    }
}

static void proj (float m[16], float a, float r, float n, float f) {
    m[5] = 1.0f / tanf (a * 0.5f);
    m[0] = m[5] / r;
    m[10] = -f / (f - n * 2.0f);
    m[11] = -2.0f * n * (f / (f - n));
    m[14] = -1.0f;

    m[1] = m[2] = m[3] = m[4] = m[6] = m[7] =
        m[8] = m[9] = m[12] = m[13] = m[15] = 0.0f;
}


static void transpose (float m[16]) {
    float t;
    t = m[1]; m[1] = m[4]; m[4] = t;
    t = m[2]; m[2] = m[8]; m[8] = t;
    t = m[3]; m[3] = m[12]; m[12] = t;
    t = m[6]; m[6] = m[9]; m[9] = t;
    t = m[7]; m[7] = m[13]; m[13] = t;
    t = m[11]; m[11] = m[14]; m[14] = t;
}

#define SCREEN_W 640
#define SCREEN_H 480

#define RAD (0.0174532925)
static void setup_view (int rx, int ry, int dist) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef (0.0, 0.0, -dist);
    glRotatef (rx, 1.0, 0.0, 0.0);
    glRotatef (ry, 0.0, 0.0, 1.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float matrix[16];
    proj (matrix, 70.0 * RAD, (float)SCREEN_W / SCREEN_H, 0.1, 1000.0);
    transpose (matrix);
    glLoadMatrixf (matrix);
}

#define GRID_W 15
#define GRID_H 15
#define n_vert (GRID_W * GRID_H)
#define n_ind ((GRID_W - 1) * (GRID_H - 1) * 2 * 3)

static void mk_grid (int *vertices, int *indices, int *colors) {
    int i, j;

    for (i = 0; i < GRID_H; i++) {
        for (j = 0; j < GRID_W; j++) {
            int index = i * GRID_W + j;
            vertices[index * 3] = j;
            vertices[index * 3 + 1] = i;
            vertices[index * 3 + 2] = (int)(3.0 * sin(i * 0.69 - 0.5));
            /* vertices[index * 3 + 2] = 0; */
            colors[index] = 0;
            if (j == 0)
                colors[index] |= MX;
            if (j == GRID_W - 1)
                colors[index] |= PX;
            if (i == 0)
                colors[index] |= MY;
            if (i == GRID_W - 1)
                colors[index] |= PY;
        }
    }
    for (i = 0; i < n_ind / 6; i++) {
        int k = i + i / (GRID_W - 1);
        indices[i * 6] = k;
        indices[i * 6 + 1] = k + 1;
        indices[i * 6 + 2] = k + GRID_W;
        indices[i * 6 + 3] = k + 1;
        indices[i * 6 + 4] = k + GRID_W + 1;
        indices[i * 6 + 5] = k + GRID_W;
    }
}

static void draw (int *vertices, int *indices, int *color, int n_indices) {
    int coul[6] = {0x00FF0000,
                   0x00FFF000,
                   0x0000FF00,
                   0x0000FFF0,
                   0x000000FF,
                   0x00F000FF};
    glPolygonMode(GL_FRONT, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < n_indices; i++) {
        int k = indices[i];
        int *vertex = &vertices[k * 3];
        int c = 0x00555555;
        if (color) {
            for (int j = 0; j < 6; j++) {
                if (color[k] & 1 << j)
                    c |= coul[j];
            }
        }
        glColor3ub (c >> 16, c >> 8 & 0xFF, c & 0xFF);
        glVertex3f (vertex[0], vertex[1], vertex[2]);
    }
    glEnd();
}


static int follow_merge (int *v_merge, int i) {
    while (v_merge[i] > -1)
        i = v_merge[i];
    return i;
}

static void reduce_mesh (int *vertices, int *indices, int *v_merge,
                         int *n_vertices, int *n_indices) {
    int i;
    int *offset = malloc (*n_vertices * sizeof *offset);
    int off = 0;

    memset (offset, 0, *n_vertices);
    for (i = 0; i < *n_vertices; i++) {
        /* offset[i] = 0; */
        while (i + off < *n_vertices && v_merge[i + off] > -1)
            off++;
        if (i + off < *n_vertices) {
            offset[i + off] = off;
            vertices[i * 3]     = vertices[(i + off) * 3];
            vertices[i * 3 + 1] = vertices[(i + off) * 3 + 1];
            vertices[i * 3 + 2] = vertices[(i + off) * 3 + 2];
        }
    }

    for (i = 0; i < *n_indices; i++) {
        int k = indices[i];
        int new = follow_merge (v_merge, k);
        indices[i] = new - offset[new];
    }
    *n_vertices -= off;
    off = 0;

    for (i = 0; i < *n_indices; i += 3) {
        while (i + off < *n_indices &&
               (indices[i + off]     == indices[i + off + 1] ||
                indices[i + off + 1] == indices[i + off + 2] ||
                indices[i + off + 2] == indices[i + off])) {
            off += 3;
        }
        if (i + off + 2 < *n_indices) {
            indices[i]     = indices[i + off];
            indices[i + 1] = indices[i + off + 1];
            indices[i + 2] = indices[i + off + 2];
        }
    }
    *n_indices -= off;

    free (offset);
}

int main (void) {
    SDL_Window *Window = NULL;
    SDL_GLContext glContext;
    const int ww = SCREEN_W, wh = SCREEN_H;

    SDL_Init (SDL_INIT_VIDEO);
    Window = SDL_CreateWindow ("LE MAO", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                               ww, wh, SDL_WINDOW_OPENGL);
    glContext = SDL_GL_CreateContext (Window);
    SDL_ShowWindow (Window);

    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glClear (GL_COLOR_BUFFER_BIT);
    glViewport (0.0, 0.0, ww, wh);
    SDL_Event ev;
    int running = 1;

    srand(1547917722);

#if 0
#define N_vertices 10
#define N_indices (N_vertices * 2 - 2)
    int n_indices = N_indices;
    int vertices[2 * N_vertices];
    int indices[N_indices];
    for (int i = 0; i < N_vertices - 1; i++) {
        indices[i * 2] = i;
        indices[i * 2 + 1] = i + 1;
    }
    for (int i = 0; i < N_vertices; i++) {
        vertices[i * 2] = i * 3 + ran (0, 2) - (N_vertices * 3 / 2);
        vertices[i * 2 + 1] = ran (-N_vertices, N_vertices);
    }
    
    int *merge = malloc (N_vertices * sizeof *merge);
    decimate (N_vertices, N_indices, vertices, indices, merge, PZ);

    /* update indices */
    for (int i = 0; i < n_indices; i++) {
        if (merge[indices[i]] > -1)
            indices[i] = merge[indices[i]];
    }
    for (int i = 0; i < n_indices; i += 2) {
        if (indices[i] == indices[i + 1]) {
            if (i < (n_indices - 1)) {
                indices[i] = indices[n_indices - 2];
                indices[i + 1] = indices[n_indices - 1];
            }
            n_indices -= 2;
            /* printf ("removed ONE line\n"); */
        }
    }
#else
    glEnable(GL_DEPTH_TEST);

    int grid_vertices[n_vert * 3];
    int grid_indices[n_ind];
    int grid_colors[n_vert];
    int grid_n_vertices = n_vert, grid_n_indices = n_ind;
    int v_merge[n_vert];
    int *v_flags = grid_colors;

    for (int i = 0; i < n_vert; i++)
        v_merge[i] = -1;
    mk_grid (grid_vertices, grid_indices, grid_colors);

    int tmp_indices[n_ind];
    int k;
    for (int i = 0; i < 6; i++) {
        char side = 1 << i;
        k = 0;
        /* capture every triangle's edge that's on side[i] */
        for (int j = 0; j < grid_n_indices; j += 3) {
            /* edge 1 */
            /* or :
               v_flags[indices[j]] & v_flags[indices[j + 1]] & side */
            if (v_flags[grid_indices[j]] & side && v_flags[grid_indices[j + 1]] & side) {
                tmp_indices[k] = grid_indices[j];
                tmp_indices[k + 1] = grid_indices[j + 1];
                k += 2;
            } else if (v_flags[grid_indices[j + 1]] & side && v_flags[grid_indices[j + 2]] & side) {
                tmp_indices[k] = grid_indices[j + 1];
                tmp_indices[k + 1] = grid_indices[j + 2];
                k += 2;
            } else if (v_flags[grid_indices[j + 2]] & side && v_flags[grid_indices[j]] & side) {
                tmp_indices[k] = grid_indices[j + 2];
                tmp_indices[k + 1] = grid_indices[j];
                k += 2;
            }
        }
        /* decimate */
        /* TODO: use v_flags to know which vertices shouldnt be touched/moved */
        decimate (grid_n_vertices, k, grid_vertices, tmp_indices, v_merge, side);
    }

    printf ("before n indices : %d\n", grid_n_indices);
    reduce_mesh (grid_vertices, grid_indices, v_merge, &grid_n_vertices, &grid_n_indices);
    printf ("after n indices : %d\n", grid_n_indices);
#endif

    int prev_x = 0, prev_y = 0, ry = 0, rx = 0, mouse_pressed = 0;
    float dist = 10.0;
    while (running){
        while (SDL_PollEvent (&ev)) {
            switch (ev.type) {
            case SDL_QUIT:
                running = 0;
                break;
            case SDL_KEYDOWN:
                switch (ev.key.keysym.sym) {
				case SDLK_ESCAPE:
					running = 0;
					break;
				default:
					break;
                }
                break;

            case SDL_MOUSEBUTTONDOWN:
                prev_x = ev.button.x;
                prev_y = ev.button.y;
                mouse_pressed = 1;
                break;

            case SDL_MOUSEBUTTONUP:
                mouse_pressed = 0;
                break;

            case SDL_MOUSEWHEEL:
                dist -= ev.wheel.y;
                break;

            case SDL_MOUSEMOTION:
                if (mouse_pressed) {
                    ry += ev.motion.x - prev_x;
                    rx += ev.motion.y - prev_y;
                    prev_x = ev.motion.x;
                    prev_y = ev.motion.y;
                }
                break;
            }
        }
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        setup_view (rx, ry, dist);

#if 1
        draw (grid_vertices, grid_indices, NULL, grid_n_indices);
#else
        glBegin (GL_LINES);
        float f = 0.5 / N_vertices;
        for (int i = 0; i < n_indices; i++) {
            glColor3f (1.0, (float)i / N_indices, f * vertices[indices[i] * 2 + 1]);
            glVertex3f(f * vertices[indices[i] * 2],
                       f * vertices[indices[i] * 2 + 1], 0.0f);
        }
        glEnd();
#endif

        SDL_GL_SwapWindow (Window);
    }

    return 0;

}
