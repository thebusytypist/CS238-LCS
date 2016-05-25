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

int Score(
    const int* prev, const int* cur,
    const char* u, const char* v,
    int i, int j,
    int step,
    int yb0, int yb1) {
    bool assigned = false;
    int opt = 0, y;

    y = j - step;
    if (y >= yb0 && y < yb1) {
        int s = cur[y] + SPScore('_', v[j]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
        }
    }

    if (prev) {
        y = j;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        y = j - step;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], v[j]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }
    } // if (prev)

    return opt;
}

void Step(
    const int* prev, int* cur,
    const char* u, const char* v,
    int i,
    int y0, int y1,
    int step,
    int yb0, int yb1) {
    for (int j = y0; j != y1; j += step) {
        cur[j] = Score(prev, cur, u, v, i, j, step, yb0, yb1);
    }
}

void Initialize(
    int* cur,
    const char* u, const char* v,
    int y0, int y1, int step) {
    cur[y0] = SCORE_SPACES;
    for (int j = y0 + 1; j != y1; j += step) {
        cur[j] = SCORE_INDEL * (abs(j - y0) + 1);
    }
}

struct Context {
    int* prev;
    int* left;
    int* right;
    int* path;

    const char* u;
    const char* v;
};

void Solve(Context* ctx, int x0, int y0, int x1, int y1) {
    if (x1 - x0 <= 1)
        return;

    int m = (x0 + x1) / 2;

    Initialize(ctx->left, ctx->u, ctx->v, y0, y1, 1);
    for (int i = x0 + 1; i <= m; ++i) {
        Step(
            ctx->left, ctx->prev, ctx->u, ctx->v,
            i,
            y0, y1,
            1,
            y0, y1);
        swap(ctx->left, ctx->prev);
    }

    Initialize(ctx->right, ctx->u, ctx->v, y1, y0 - 1, -1);
    for (int i = x1 - 1; i > m; --i) {
        Step(
            ctx->right, ctx->prev, ctx->u, ctx->v,
            i,
            y1, y0 - 1,
            -1,
            y0, y1 + 1);
        swap(ctx->right, ctx->prev);
    }

    int k = y0, l = ctx->left[y0] + ctx->right[y0 + 1];
    for (int j = y0 + 1; j < y1; ++j) {
        int n = ctx->left[j] + ctx->right[j + 1];
        if (n > l) {
            l = n;
            k = j;
        }
    }

    ctx->path[m] = k;

    Solve(ctx, x0, y0, m, k);
    Solve(ctx, m, k, x1, y1);
}

int main(int argc, char* argv[]) {
    char str0[MAXSTRLEN], str1[MAXSTRLEN];

    scanf("%s%s", str0 + 1, str1 + 1);
    str0[0] = str1[0] = '_';
    const char* u = str0;
    const char* v = str1;
    int lu = strlen(u), lv = strlen(v);

    if (lu > lv) {
        swap(lu, lv);
        swap(u, v);
    }

    int buffer[4 * MAXSTRLEN];

    Context ctx;

    ctx.prev = buffer;
    ctx.left = ctx.prev + MAXSTRLEN;
    ctx.right = ctx.left + MAXSTRLEN;
    ctx.path = ctx.right + MAXSTRLEN;
    ctx.u = u;
    ctx.v = v;

    Solve(&ctx, 0, 0, lu, lv);

    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.path[i]);
    }
    printf("\n");

    //--------------------------------------------------------------------------
    int total = 0, length = 0, matches = 0;
    for (int i = 1, b = ctx.path[1]; i < lu; ++i) {
        while (b < ctx.path[i]) {
            total += SPScore('_', v[b]);
            ++length;
            ++b;
        }
        total += SPScore(u[i], v[b]);
        matches += (u[i] == v[b]);
        ++length;
        b = ctx.path[i] + 1;
    }
    if (ctx.path[lu - 1] + 1 < lv) {
        int t = lv - (ctx.path[lu - 1] + 1);
        total += SCORE_INDEL * t;
        length += t;
    }

    printf("score: %d\nlength: %d\nmatches: %d\n", total, length, matches);

    //--------------------------------------------------------------------------
    Initialize(ctx.left, ctx.u, ctx.v, 0, lv, 1);
    for (int i = 1; i < lu; ++i) {
        Step(
            ctx.left, ctx.prev, ctx.u, ctx.v,
            i,
            0, lv,
            1,
            0, lv);
        swap(ctx.left, ctx.prev);
    }
    printf("reference score: %d\n", ctx.left[lv - 1]);

    //--------------------------------------------------------------------------
    for (int i = 1, b = ctx.path[1]; i < lu; ++i) {
        while (b < ctx.path[i]) {
            printf("_ ");
            ++b;
        }
        b = ctx.path[i] + 1;

        printf("%c ", u[i]);
    }
    for (int i = ctx.path[lu - 1] + 1; i < lv; ++i) {
        printf("_ ");
    }
    printf("\n");

    for (int i = 1, b = ctx.path[1]; i < lu; ++i) {
        while (b <= ctx.path[i]) {
            printf("%c ", v[b]);
            ++b;
        }
    }
    for (int b = ctx.path[lu - 1] + 1; b < lv; ++b) {
        printf("%c ", v[b]);
    }
    printf("\n");

    return 0;
}
