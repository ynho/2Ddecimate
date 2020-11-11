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

#define EPSILON 0.001

#define TRUE 1
#define FALSE 0

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
    float area, distance;
    int a, b;
};

struct mesh {
    int *vertices;
    int *indices;
    int *v_flags;
    int n_vertices;
    int n_indices;
    int v_edgestart;
    int i_edgestart;
};

static void init_mesh (struct mesh *mesh) {
    mesh->vertices = NULL;
    mesh->indices = NULL;
    mesh->v_flags = NULL;
    mesh->n_vertices = 0;
    mesh->n_indices = 0;
    mesh->v_edgestart = 0;
    mesh->i_edgestart = 0;
}

#if 0
static int ran (int a, int b) {
    float r = (float)rand () / RAND_MAX;
    r *= (b - a);
    r += a;
    return r + 0.5;
}
#endif

static float distance (int x[2], int y[2]) {
    int a = x[0] - y[0];
    int b = x[1] - y[1];
    return sqrt (a * a + b * b);
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

static void coordinates_from_side (int side, int *x, int *y) {
    if (side & (MX | PX)) {
        *x = 1, *y = 2;
    } else if (side & (MY | PY)) {
        *x = 0, *y = 2;
    } else {
        *x = 0, *y = 1;
    }
}

static void compute_area_and_distance (struct tri *tri, int self, int *vertices, int side) {
    int x, y;
    coordinates_from_side (side, &x, &y);
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
        tri->distance = distance (b, c);
    }
}

static void compute_areas (struct tri *origin, struct tri **tris, int n_tris,
                           int *vertices, int side) {
    for (int i = 0; i < n_tris; i++) {
        int index = (int)(tris[i] - origin);
        compute_area_and_distance (tris[i], index, vertices, side);
    }
}

#if 0
/* wtf is this? */
static int secondary_condition(int *vertices, int a, int b) {
    if (vertices[a * 3] == vertices[b * 3]) {
        if (vertices[a * 3 + 1] == vertices[b * 3 + 1])
            return vertices[a * 3 + 2] < vertices[b * 3 + 2];
        else
            return vertices[a * 3 + 1] < vertices[b * 3 + 1];
    } else
        return vertices[a * 3] < vertices[b * 3];
}
#endif

static void sort_areas (struct tri **tris, int *vertices, int n_tris) {
    struct tri *x = NULL;
    for (int i = 1; i < n_tris; i++) {
        int j = i;
        while (j >= 1 &&
               ((fabs (tris[j]->area - tris[j - 1]->area) < EPSILON
                 && tris[j]->distance < tris[j - 1]->distance)
                ||
                (tris[j]->area < tris[j - 1]->area
                 && fabs (tris[j]->area - tris[j - 1]->area) >= EPSILON
                /* wtf is this? */
                /* || */
                /* (tris[j]->area == tris[j - 1]->area && */
                /*  secondary_condition(vertices, j, j - 1)) */
                    ))) {
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

void decimate (int n_vertices, int n_indices, int *vertices, int *indices,
               int *merge, int side, struct tri *tris, struct tri **sorted) {
    int i;
    size_t n_sorted = n_indices / 2; /* divided by 2 because 1 edge has 2 vertices */
    int offset = 0;

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

    int reduces = n_indices / 4;

    for (i = 0; i < reduces; i++) {
        sort_areas (&sorted[offset], vertices, n_sorted - offset);

        offset = 0;
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

        compute_area_and_distance (&tris[a], a, vertices, side);
        compute_area_and_distance (&tris[b], b, vertices, side);

        n_sorted--;
        shift_left (sorted, offset, n_sorted);
    }
}

static int follow_merge (int *v_merge, int i) {
    while (v_merge[i] > -1)
        i = v_merge[i];
    return i;
}

static void reduce_vertices (struct mesh *mesh, int *v_merge, int *offset) {
    int i, off = 0, n_edge_vertices = mesh->n_vertices - mesh->v_edgestart;
    memset (offset, 0, n_edge_vertices * sizeof *offset);
    for (i = 0; i < n_edge_vertices; i++) {
        /* offset[i] = 0; */
        while (i + off < n_edge_vertices && v_merge[i + off] > -1)
            off++;
        if (i + off < n_edge_vertices) {
            int j = i + mesh->v_edgestart;
            offset[i + off] = off;
            mesh->vertices[j * 3]     = mesh->vertices[(j + off) * 3];
            mesh->vertices[j * 3 + 1] = mesh->vertices[(j + off) * 3 + 1];
            mesh->vertices[j * 3 + 2] = mesh->vertices[(j + off) * 3 + 2];
            mesh->v_flags[j] = mesh->v_flags[j + off];
        }
    }
    mesh->n_vertices -= off;
}

static void reduce_mesh (struct mesh *mesh, int *v_merge, int *offset, int reduce_indices) {
    int i;

    reduce_vertices (mesh, v_merge, offset);

    for (i = mesh->i_edgestart; i < mesh->n_indices; i++) {
        if (mesh->indices[i] >= mesh->v_edgestart) {
            int new = follow_merge (v_merge, mesh->indices[i] - mesh->v_edgestart);
            mesh->indices[i] = new - offset[new] + mesh->v_edgestart;
        }
    }

    if (reduce_indices) {
        int off = 0;
        for (i = mesh->i_edgestart; i < mesh->n_indices; i += 3) {
            while (i + off < mesh->n_indices &&
                   (mesh->indices[i + off]     == mesh->indices[i + off + 1] ||
                    mesh->indices[i + off + 1] == mesh->indices[i + off + 2] ||
                    mesh->indices[i + off + 2] == mesh->indices[i + off])) {
                off += 3;
            }
            if (i + off + 2 < mesh->n_indices) {
                mesh->indices[i]     = mesh->indices[i + off];
                mesh->indices[i + 1] = mesh->indices[i + off + 1];
                mesh->indices[i + 2] = mesh->indices[i + off + 2];
            }
        }
        mesh->n_indices -= off;
    }
}


void full_decimate_edges (int edges, struct mesh *mesh) {
    int i, k;
    int *indices = mesh->indices, *v_flags = mesh->v_flags;
    /* int voff = mesh->v_edgestart; */
    /* int ioff = mesh->v_edgestart; */
    size_t i2, i3, i4;
    size_t n_edge_vertices = mesh->n_vertices - mesh->v_edgestart;
    size_t n_edge_indices = mesh->n_indices - mesh->i_edgestart;
    i2 =      n_edge_vertices * sizeof (int);
    i3 = i2 + n_edge_indices  * sizeof (int);
    i4 = i3 + n_edge_vertices * sizeof (struct tri);
    size_t total = i4 + n_edge_indices * sizeof (struct tri*);
    /* TODO: optimize memory allocation by reusing memory */
    void *memory = malloc (total);
    int *v_merge = memory;
    int *tmp_indices = memory + i2;
    struct tri *tris = memory + i3;
    struct tri **sorted = memory + i4;

    for (i = 0; i < n_edge_vertices; i++)
        v_merge[i] = -1;

    for (i = 0; i < 6; i++) {
        char side = 1 << i;
        if (!(side & edges))
            continue;
        k = 0;
        /* capture every triangle's edge that's on side[i] */
        for (int j = mesh->i_edgestart; j < mesh->n_indices; j += 3) {
            /* or :
               v_flags[indices[j]] & v_flags[indices[j + 1]] & side */
            /* TODO: v_flags should only apply from v_edgestart */
            if        (v_flags[indices[j]]     & side && v_flags[indices[j + 1]] & side) {
                tmp_indices[k]     = indices[j]     - mesh->v_edgestart;
                tmp_indices[k + 1] = indices[j + 1] - mesh->v_edgestart;
                k += 2;
            } else if (v_flags[indices[j + 1]] & side && v_flags[indices[j + 2]] & side) {
                tmp_indices[k]     = indices[j + 1] - mesh->v_edgestart;
                tmp_indices[k + 1] = indices[j + 2] - mesh->v_edgestart;
                k += 2;
            } else if (v_flags[indices[j + 2]] & side && v_flags[indices[j]]     & side) {
                tmp_indices[k]     = indices[j + 2] - mesh->v_edgestart;
                tmp_indices[k + 1] = indices[j]     - mesh->v_edgestart;
                k += 2;
            }
        }
        /* decimate */
        /* TODO: use v_flags to know which vertices shouldnt be touched/moved */
        decimate (mesh->n_vertices - mesh->v_edgestart, k, &mesh->vertices[mesh->v_edgestart * 3],
                  tmp_indices, v_merge, side, tris, sorted);
    }

    /* TODO: move memory allocation */
    int *offset = malloc (mesh->n_vertices * sizeof *offset);
    reduce_mesh (mesh, v_merge, offset, TRUE);
    free (offset);

    free (memory);
}


static void concatenate_meshes (struct mesh *a, struct mesh *b, struct mesh *result) {
    int i;

    result->n_vertices = a->n_vertices + b->n_vertices;
    result->vertices = malloc (result->n_vertices * sizeof *result->vertices * 3);
    result->v_flags = malloc (result->n_vertices * sizeof *result->v_flags);
    result->n_indices = a->n_indices + b->n_indices;
    result->indices = malloc (sizeof *result->indices * result->n_indices);

    memcpy(result->vertices, a->vertices, a->n_vertices * 3 * sizeof *result->vertices);
    memcpy(&result->vertices[a->n_vertices * 3], b->vertices, b->n_vertices * 3 * sizeof *result->vertices);
    memcpy(result->v_flags, a->v_flags, a->n_vertices * sizeof *result->v_flags);
    memcpy(&result->v_flags[a->n_vertices], b->v_flags, b->n_vertices * sizeof *result->v_flags);
    memcpy(result->indices, a->indices, a->n_indices * sizeof *result->indices);
    memcpy(&result->indices[a->n_indices], b->indices, b->n_indices * sizeof *result->indices);
    for (i = 0; i < b->n_indices; i++)
        result->indices[i + a->n_indices] += a->n_vertices;
}

static int compare_vertices (int side, int *v1, int *v2) {
    int x, y;
    coordinates_from_side (side, &x, &y);
    return v1[x] == v2[x] && v1[y] == v2[y];
}

void merge (int side, struct mesh *a, struct mesh *b, struct mesh *result) {
    int i, j;
    int *offset = NULL;
    int *v_merge = NULL;

    concatenate_meshes (a, b, result);

    v_merge = malloc (result->n_vertices * sizeof *v_merge);
    offset = malloc (result->n_vertices * sizeof *offset);
    for (i = 0; i < result->n_vertices; i++)
        v_merge[i] = -1;

    for (i = a->n_vertices; i < result->n_vertices; i++) {
        if (result->v_flags[i] & (side << 1)) {
            /* find matching vertex */
            /* NOTE: make sure not to overlap with B's vertices, or there would be self-merging */
            for (j = 0; j < a->n_vertices; j++) {
                if (result->v_flags[j] & side &&
                    compare_vertices (side, &result->vertices[j * 3], &result->vertices[i * 3])) {
                    v_merge[i] = j;
                    result->v_flags[j] ^= side;
                }
            }
        }
    }

    reduce_mesh (result, v_merge, offset, FALSE);

    free (v_merge);
    free (offset);
}

int face_mask_from_cube_id (int id) {
    return MX << (id & 1) | MY << (id & 2 ? 1 : 0) | MZ << (id & 4 ? 1 : 0);
}

void merge_cubes (struct mesh m[8], struct mesh *result) {
    int i;
    struct mesh mesh[6];

    /* decimate _outside_ faces to match with lower-level LODs */
    for (i = 0; i < 8; i++)
        full_decimate_edges (face_mask_from_cube_id (i), &m[i]);

    for (i = 0; i < 4; i++)
        merge (MX, &m[i * 2 + 1], &m[i * 2], &mesh[i]);
    merge (MY, &mesh[1], &mesh[0], &mesh[4]);
    merge (MY, &mesh[3], &mesh[2], &mesh[5]);
    merge (MZ, &mesh[5], &mesh[4], result);
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

#define SCREEN_W 1600
#define SCREEN_H 900

#define RAD (0.0174532925)
static void setup_view (int rx, int ry, int dist) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef (0.0, 0.0, -dist);
    glRotatef (rx, 1.0, 0.0, 0.0);
    glRotatef (ry, 0.0, 0.0, 1.0);
    glTranslatef (-10.0, -10.0, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float matrix[16];
    proj (matrix, 70.0 * RAD, (float)SCREEN_W / SCREEN_H, 0.1, 1000.0);
    transpose (matrix);
    glLoadMatrixf (matrix);
}

#define GRID_W 15
#define GRID_H 16
#define n_vert (GRID_W * GRID_H)
#define n_ind ((GRID_W - 1) * (GRID_H - 1) * 2 * 3)

static void mk_grid (int x, int y, int z, struct mesh *mesh) {
    int i, j;

    for (i = 0; i < GRID_H; i++) {
        for (j = 0; j < GRID_W; j++) {
            int index = i * GRID_W + j;
            mesh->vertices[index * 3] = j + x;
            mesh->vertices[index * 3 + 1] = i + y;
            mesh->vertices[index * 3 + 2] = (int)(3.0 * sin(i * 0.69 - 0.5)) + z;
            /* vertices[index * 3 + 2] = 0; */
            mesh->v_flags[index] = 0;
            if (j == 0)
                mesh->v_flags[index] |= MX;
            if (j == GRID_W - 1)
                mesh->v_flags[index] |= PX;
            if (i == 0)
                mesh->v_flags[index] |= MY;
            if (i == GRID_H - 1)
                mesh->v_flags[index] |= PY;
        }
    }
    for (i = 0; i < n_ind / 6; i++) {
        int k = i + i / (GRID_W - 1);
        mesh->indices[i * 6] = k;
        mesh->indices[i * 6 + 1] = k + 1;
        mesh->indices[i * 6 + 2] = k + GRID_W;
        mesh->indices[i * 6 + 3] = k + 1;
        mesh->indices[i * 6 + 4] = k + GRID_W + 1;
        mesh->indices[i * 6 + 5] = k + GRID_W;
    }
}

static inline void switch_vertices (int *a, int *b) {
    int c;
    c = a[0], a[0] = b[0], b[0] = c;
    c = a[1], a[1] = b[1], b[1] = c;
    c = a[2], a[2] = b[2], b[2] = c;
}

/* more accurately: switch triangles */
static inline void switch_indices (int *a, int *b) {
    return switch_vertices (a, b);
}

static inline void switch_v_flags (int *a, int *b) {
    int c;
    c = *a, *a = *b, *b = c;
}

static inline int is_triangle_on_edge (int *indices, int *v_flags) {
    return v_flags[indices[0]] || v_flags[indices[1]] || v_flags[indices[2]];
}

static size_t partition_vertices (int num, int *vertices, int *v_flags, int *map) {
    int i, end = num - 1;
    for (i = 0; i < num; i++)
        map[i] = -1;
    for (i = 0; i < end; i++) {
        if (v_flags[i]) {
            switch_vertices (&vertices[i * 3], &vertices[end * 3]);
            switch_v_flags (&v_flags[i], &v_flags[end]);
            map[end] = i;
            if (end < num - 1 && map[end + 1] == i)
                map[end + 1] = end;
            else
                map[i] = end;
            i--;
            end--;
        }
    }
    return v_flags[i] ? i : i + 1;
}

static size_t partition_indices (int num, int *indices, int *v_flags) {
    int i, end = num - 1;
    for (i = 0; i < end; i++) {
        if (is_triangle_on_edge(&indices[i * 3], v_flags)) {
            switch_indices (&indices[i * 3], &indices[end * 3]);
            i--;
            end--;
        }
    }
    return is_triangle_on_edge(&indices[i * 3], v_flags) ? i : i + 1;
}

static void separate_edges (struct mesh *mesh) {
    int i;
    int *map = malloc (mesh->n_vertices * sizeof *map);

    mesh->v_edgestart = partition_vertices (mesh->n_vertices, mesh->vertices, mesh->v_flags, map);
    for (i = 0; i < mesh->n_indices; i++) {
        if (map[mesh->indices[i]] > -1)
            mesh->indices[i] = map[mesh->indices[i]];
    }
    mesh->i_edgestart = partition_indices (mesh->n_indices / 3, mesh->indices, mesh->v_flags);
    mesh->i_edgestart *= 3;
}

static void draw (struct mesh *mesh) {
    int coul[6] = {0x00550000,
                   0x00555000,
                   0x00005500,
                   0x00005550,
                   0x00000055,
                   0x00500055};
    glPolygonMode(GL_FRONT, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < mesh->n_indices; i++) {
        int k = mesh->indices[i];
        int *vertex = &mesh->vertices[k * 3];
        int c = 0x00222222;
        if (mesh->v_flags) {
            for (int j = 0; j < 6; j++) {
                if (mesh->v_flags[k] & (1 << j)) {
                    c += coul[j];
                }
            }
        }
        glColor3ub (c >> 16, c >> 8 & 0xFF, c & 0xFF);
        glVertex3f (vertex[0], vertex[1], vertex[2]);
    }
    glEnd();
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

    glEnable(GL_DEPTH_TEST);

    int i;
    int grid_vertices[8][n_vert * 3];
    int grid_indices[8][n_ind];
    int grid_v_flags[8][n_vert];
    struct mesh mesh[9];

    for (i = 0; i < 9; i++)
        init_mesh (&mesh[i]);

    for (i = 0; i < 8; i++) {
        mesh[i].vertices = grid_vertices[i];
        mesh[i].indices = grid_indices[i];
        mesh[i].v_flags = grid_v_flags[i];
        mesh[i].n_vertices = n_vert;
        mesh[i].n_indices = n_ind;
    }

    mk_grid (0, 0, 0, &mesh[0]);
    /* mk_grid (GRID_W, 0, 0, &mesh[1]); */
    /* mk_grid (0, GRID_H, 0, &mesh[2]); */
    /* mk_grid (GRID_W, GRID_H, 0, &mesh[3]); */
    /* merge_cubes (mesh, &mesh[8]); */

    separate_edges (&mesh[0]);
    full_decimate_edges (MX | PX | MY | PY, &mesh[0]);

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

        draw (&mesh[0]);

        SDL_GL_SwapWindow (Window);
    }

    return 0;

}
