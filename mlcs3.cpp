#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::max;
using std::abs;
using std::swap;

const int MAXSTRLEN = 128;

const int SCORE_MATCH = 5;
const int SCORE_MISMATCH = -4;
const int SCORE_INDEL = -8;
const int SCORE_SPACES = 0;

struct Context {
    int* prev;
    int* left;
    int* right;

    int* pathv;
    int* pathw;

    const char* u;
    const char* v;
    const char* w;

    int dimv, dimw;
};

int SPScore(char a, char b) {
    if (a == b && a != '_') {
        return SCORE_MATCH;
    }
    else if (a != b && a != '_' && b != '_') {
        return SCORE_MISMATCH;
    }
    else if (a != '_' || b != '_') {
        return SCORE_INDEL;
    }
    else
        return SCORE_SPACES;
}

int SPScore(char a, char b, char c) {
    return SPScore(a, b) + SPScore(b, c) + SPScore(c, a);
}

int Score(
    const Context* ctx,
    const int* prev, const int* cur,
    int i, int j, int k,
    int step,
    int yb0, int yb1,
    int zb0, int zb1) {
    int dimy = ctx->dimv, dimz = ctx->dimw;
    bool assigned = false;
    int opt = 0, y, z;

    y = j - step;
    z = k;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', ctx->v[j], '_');
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
        }
    }

    y = j;
    z = k - step;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', '_', ctx->w[k]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
        }
    }

    y = j - step;
    z = k - step;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', ctx->v[j], ctx->w[k]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
        }
    }

    if (prev) {
        y = j;
        z = k;
        {
            int s = prev[y * dimz + z] + SPScore(ctx->u[i], '_', '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        y = j - step;
        z = k;
        if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
            int s = prev[y * dimz + z] + SPScore(ctx->u[i], ctx->v[j], '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        y = j;
        z = k - step;
        if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
            int s = prev[y * dimz + z] + SPScore(ctx->u[i], '_', ctx->w[k]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        y = j - step;
        z = k - step;
        if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
            int s = prev[y * dimz + z]
                + SPScore(ctx->u[i], ctx->v[j], ctx->w[k]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }
    } // if (prev)

    return opt;
}

void Step(
    const Context* ctx,
    const int* prev, int* cur,
    int i,
    int y0, int z0,
    int y1, int z1,
    int step,
    int yb0, int yb1,
    int zb0, int zb1) {
    int dimz = ctx->dimw;
    for (int j = y0; j != y1; j += step) {
        for (int k = z0; k != z1; k += step) {
            cur[j * dimz + k] =
                Score(ctx, prev, cur, i, j, k, step, yb0, yb1, zb0, zb1);
        }
    }
}

void Initialize(
    const Context* ctx,
    int* cur,
    int y0, int z0,
    int y1, int z1,
    int step,
    int yb0, int yb1,
    int zb0, int zb1) {
    int dimy = ctx->dimv, dimz = ctx->dimw;
    cur[y0 * dimz + z0] = SCORE_SPACES;
    for (int k = z0 + 1; k != z1; k += step) {
        cur[y0 * dimz + k] = 2 * SCORE_INDEL * abs(k - z0);
    }

    for (int j = y0 + 1; j != y1; j += step) {
        for (int k = z0; k != z1; k += step) {
            cur[j * dimz + k] =
                // u[0] is '_'.
                Score(ctx, nullptr, cur, 0, j, k, step, yb0, yb1, zb0, zb1);
        }
    }
}

void Solve(Context* ctx, int x0, int y0, int z0, int x1, int y1, int z1) {
    if (x1 - x0 <= 1)
        return;

    int m = (x0 + x1) / 2;

    Initialize(ctx, ctx->left, y0, z0, y1, z1, 1, y0, y1, z0, z1);
    for (int i = x0 + 1; i <= m; ++i) {
        Step(
            ctx, ctx->left, ctx->prev,
            i,
            y0, z0,
            y1, z1,
            1,
            y0, y1,
            z0, z1);
        swap(ctx->left, ctx->prev);
    }

    Initialize(ctx, ctx->right,
        y1, z1, y0, z0,
        -1,
        y0 + 1, y1 + 1,
        z0 + 1, z1 + 1);
    for (int i = x1 - 1; i > m; --i) {
        Step(
            ctx, ctx->right, ctx->prev,
            i,
            y1, z1,
            y0, z0,
            -1,
            y0 + 1, y1 + 1,
            z0 + 1, z1 + 1);
        swap(ctx->right, ctx->prev);
    }

    int dimz = ctx->dimw;
    int my = y0, mz = z0;
    int l = ctx->left[my * dimz + mz]
        + ctx->right[(my + 1) * dimz + mz + 1];
    for (int j = y0; j < y1; ++j) {
        for (int k = z0; k < z1; ++k) {
            int n = ctx->left[j * dimz + k]
                + ctx->right[(j + 1) * dimz + k + 1];
            if (n > l) {
                l = n;
                my = j;
                mz = k;
            }
        }
    }

    ctx->pathv[m] = my;
    ctx->pathw[m] = mz;

    Solve(ctx, x0, y0, z0, m, my, mz);
    Solve(ctx, m, my, mz, x1, y1, z1);
}

int main() {
    char str0[MAXSTRLEN], str1[MAXSTRLEN], str2[MAXSTRLEN];
    scanf("%s%s%s", str0 + 1, str1 + 1, str2 + 1);
    str0[0] = str1[0] = str2[0] = '_';

    const char* u = str0;
    const char* v = str1;
    const char* w = str2;

    int lu = strlen(u), lv = strlen(v), lw = strlen(w);

    if (lv < lu) {
        swap(u, v);
        swap(lu, lv);
    }
    if (lw < lu) {
        swap(u, w);
        swap(lu, lw);
    }
    if (lv > lw) {
        swap(v, w);
        swap(lv, lw);
    }

    printf("%s %s %s\n", u + 1, v + 1, w + 1);

    int buffer[(MAXSTRLEN * 3 + 2) * MAXSTRLEN];

    Context ctx;
    ctx.prev = buffer;
    ctx.left = ctx.prev + MAXSTRLEN * MAXSTRLEN;
    ctx.right = ctx.left + MAXSTRLEN * MAXSTRLEN;
    ctx.pathv = ctx.right + MAXSTRLEN * MAXSTRLEN;
    ctx.pathw = ctx.pathv + MAXSTRLEN;

    ctx.u = u;
    ctx.v = v;
    ctx.w = w;

    ctx.dimv = lv;
    ctx.dimw = lw;

    //--------------------------------------------------------------------------

    Solve(&ctx, 0, 0, 0, lu, lv, lw);

    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.pathv[i]);
    }
    printf("\n");
    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.pathw[i]);
    }
    printf("\n");

    //--------------------------------------------------------------------------
    Initialize(&ctx, ctx.left,
        0, 0, lv, lw,
        1,
        0, lv,
        0, lw);
    for (int i = 1; i < lu; ++i) {
        Step(&ctx, ctx.left, ctx.prev,
            i,
            0, 0,
            lv, lw,
            1,
            0, lv,
            0, lw);
        swap(ctx.left, ctx.prev);
    }

    printf("reference score: %d\n", ctx.left[(lv - 1) * lw + lw - 1]);

    //--------------------------------------------------------------------------


    return 0;
}