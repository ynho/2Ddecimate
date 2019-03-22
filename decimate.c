#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>



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

static void compute_area (struct tri *tri, int self, int *vertices) {
    if (tri->a < 0 || tri->b < 0)
        tri->area = -1.0;
    else {
        int a[2], b[2], c[2];
        a[0] = vertices[self * 2];
        a[1] = vertices[self * 2 + 1];
        b[0] = vertices[tri->a * 2];
        b[1] = vertices[tri->a * 2 + 1];
        c[0] = vertices[tri->b * 2];
        c[1] = vertices[tri->b * 2 + 1];
        tri->area = heron (a, b, c);
    }
}

static void compute_areas (struct tri *tris, int n_tris, int *vertices) {
    for (int i = 0; i < n_tris; i++)
        compute_area (&tris[i], i, vertices);
}

static void sort_areas (struct tri **tris, int n_tris) {
    struct tri *x = NULL;
    for (int i = 1; i < n_tris; i++) {
        int j = i;
        while (tris[j]->area < tris[j - 1]->area) {
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

int* decimate (int n_vertices, int n_indices, int *vertices, int *indices, int *merge) {
    int i;
    /* no. */
    /* int *merge = malloc (n_indices * sizeof *merge); */
    /* yes.. ? */
    /* int *merge = malloc (n_vertices * sizeof *merge); */

    /* so here, n_vertices is actually correct but since we are getting
       all the vertices from the whole cube, hacking in the (tweaked) indices
       number is better (smaller number) */
    /* struct tri *tris = malloc (n_vertices * sizeof *tris); */
    struct tri *tris = malloc (n_vertices * sizeof *tris);
    struct tri **sorted = malloc (n_vertices * sizeof *tris);

    for (i = 0; i < n_vertices; i++)
        merge[i] = -1;
    for (i = 0; i < n_vertices; i++) {
        tris[i].a = tris[i].b = -1;
        sorted[i] = &tris[i];
    }
    for (i = 0; i < n_indices / 2; i++) {
        int a = indices[i * 2], b = indices[i * 2 + 1];
        v_connect (tris, a, b);
        v_connect (tris, b, a);
    }

    /* compute all areas */
    compute_areas (tris, n_vertices, vertices);

    /* sort all areas */
    sort_areas (sorted, n_vertices);

    int reduces = n_vertices / 2;

    for (i = 0; i < reduces; i++) {
        /* skip borders */
        int offset = 0;
        while (sorted[offset]->area < 0.0)
            offset++;

        int origin = (int)(sorted[offset] - tris);

        /* reduce least significant vertex */
        int a = tris[origin].a, b = tris[origin].b;

        merge[origin] = a;
        /* tris[origin].area = -1.0; */

        reconnect (&tris[a], origin, b);
        reconnect (&tris[b], origin, a);

        compute_area (&tris[a], a, vertices);
        compute_area (&tris[b], b, vertices);

        n_vertices--;
        shift_left (sorted, offset, n_vertices);
        sort_areas (&sorted[offset], n_vertices - offset);
    }

    return merge;
}
#if 0
int* decimate (int n_vertices, int n_indices, int *vertices, int *indices, int *merge) {
    int i;
    /* no. */
    /* int *merge = malloc (n_indices * sizeof *merge); */
    /* yes.. ? */
    /* int *merge = malloc (n_vertices * sizeof *merge); */

    /* so here, n_vertices is actually correct but since we are getting
       all the vertices from the whole cube, hacking in the (tweaked) indices
       number is better (smaller number) */
    /* struct tri *tris = malloc (n_vertices * sizeof *tris); */
    struct tri *tris = malloc (n_vertices * sizeof *tris);
    struct tri **sorted = malloc (n_vertices * sizeof *tris);

    for (i = 0; i < n_indices; i++)
        merge[i] = -1;
    for (i = 0; i < n_vertices; i++) {
        tris[i].a = tris[i].b = -1;
        sorted[i] = &tris[i];
    }
    for (i = 0; i < n_indices / 2; i++) {
        int a = indices[i * 2], b = indices[i * 2 + 1];
        v_connect (tris, a, b);
        v_connect (tris, b, a);
    }

    /* compute all areas */
    compute_areas (tris, n_vertices, vertices);

    /* sort all areas */
    sort_areas (sorted, n_vertices);

    int reduces = n_vertices / 2;

    for (i = 0; i < reduces; i++) {
        /* skip borders */
        int offset = 0;
        while (sorted[offset]->area < 0.0)
            offset++;

        int origin = (int)(sorted[offset] - tris);

        /* reduce least significant vertex */
        int a = tris[origin].a, b = tris[origin].b;

        merge[origin] = a;
        /* tris[origin].area = -1.0; */

        reconnect (&tris[a], origin, b);
        reconnect (&tris[b], origin, a);

        compute_area (&tris[a], a, vertices);
        compute_area (&tris[b], b, vertices);

        n_vertices--;
        shift_left (sorted, offset, n_vertices);
        sort_areas (&sorted[offset], n_vertices - offset);
    }

    return merge;
}

#define SIDE_1 1
#define SIDE_2 2
#define SIDE_3 4
#define SIDE_4 8
#define SIDE_5 16
#define SIDE_6 32

int decimate_edges (float *vertices, char *v_flags,
                    int *indices, int n_indices, int n_vertice) {
    int s, i, j, k;
    int *tmp_indices = malloc (n_indices * sizeof *tmp_indices);
    int *v_merge = malloc (n_vertices * sizeof *v_merge);

    for (i = 0; i < n_vertices; i++)
        v_merge[i] = -1;

    char sides[6] = {SIDE_1, SIDE_2, SIDE_3, SIDE_4, SIDE_5, SIDE_6};
    for (i = 0; i < 6; i++) {
        char side = sides[i];
        k = 0;
        /* capture every triangle's edge that's on side[i] */
        for (j = 0; j < n_indices; j += 3) {
            /* edge 1 */
            if (v_flags[indices[j]] & side && v_flags[indices[j + 1]] & side) {
                tmp_indices[k] = indices[j];
                tmp_indices[k + 1] = indices[j + 1];
                k += 2;
            } else if (v_flags[indices[j + 1]] & side && v_flags[indices[j + 2]] & side) {
                tmp_indices[k + 1] = indices[j + 1];
                tmp_indices[k + 2] = indices[j + 2];
                k += 2;
            } else if (v_flags[indices[j + 1]] & side && v_flags[indices[j]] & side) {
                tmp_indices[k + 1] = indices[j + 1];
                tmp_indices[k] = indices[j];
                k += 2;
            }
        }
        /* decimate */
        /* TODO: use v_flags to know which vertices shouldnt be touched/moved */
        decimate (n_vertices, k, vertices, tmp_indices, v_merge);
    }

    /* now reduce according to v_merge */
    int *offset = malloc (n_vertices * sizeof *moved);
    int off = 0;
    for (i = 0; i < n_vertices; i++) {
        while (v_merge[i + off] > -1) {
            off++;
            if (i + off >= n_vertices)
                goto end;
        }
        offset[i] = off;
        vertices[i * 3]     = vertices[(i + off) * 3];
        vertices[i * 3 + 1] = vertices[(i + off) * 3 + 1];
        vertices[i * 3 + 2] = vertices[(i + off) * 3 + 2];
    }
end:
    for (i = 0; i < n_indices; i++) {
        int new = v_merge[i] > -1 ? v_merge[i] : indices[i];
        indices[i] = new - offset[i];
    }
    n_vertices -= off;
    off = 0;
    for (i = 0; i < n_indices; i += 3) {
        j = i + off;
        while (indices[i + off]     == indices[i + off + 1] ||
               indices[i + off + 1] == indices[i + off + 2] ||
               indices[i + off + 2] == indices[i + off]) {
            off += 3;
            if (i + off + 2 >= n_indices)
                goto finito;
        }
        indices[i]     = indices[i + off];
        indices[i + 1] = indices[i + off + 1];
        indices[i + 2] = indices[i + off + 2];
    }
finito:
    n_indices -= off;

    /* tadaa */
}
#endif


int main (void) {
    SDL_Window *Window = NULL;
    SDL_GLContext glContext;
    const int ww = 640, wh = 480;

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

#define N_vertices 10
#define N_indices (N_vertices * 2 - 2)
    int n_indices = N_indices;
    int vertices[2 * N_vertices];
    int indices[N_indices];
    for (int i = 0; i < N_vertices - 1; i++) {
        indices[i * 2] = i;
        indices[i * 2 + 1] = i + 1;
    }
    srand(1547917722);
    for (int i = 0; i < N_vertices; i++) {
        vertices[i * 2] = i * 3 + ran (0, 2) - (N_vertices * 3 / 2);
        vertices[i * 2 + 1] = ran (-N_vertices, N_vertices);
    }
    
    int *merge = malloc (N_vertices * sizeof *merge);
    decimate (N_vertices, N_indices, vertices, indices, merge);

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

    while (running){

        while (SDL_PollEvent (&ev)) {
            if (ev.type == SDL_QUIT)
                running = 0;
            if (ev.type == SDL_KEYDOWN ){
                switch (ev.key.keysym.sym ){
				case SDLK_ESCAPE:
					running = 0;
					break;
				default:
					break;
                }
            }
        }
        glClear (GL_COLOR_BUFFER_BIT);

        glBegin (GL_LINES);
        float f = 0.5 / N_vertices;
        for (int i = 0; i < n_indices; i++) {
            glColor3f (1.0, (float)i / N_indices, f * vertices[indices[i] * 2 + 1]);
            glVertex3f(f * vertices[indices[i] * 2],
                       f * vertices[indices[i] * 2 + 1], 0.0f);
        }
        glEnd();

        SDL_GL_SwapWindow (Window);
    }

    return 0;

}
