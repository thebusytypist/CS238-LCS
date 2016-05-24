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

int SPScore(char a, char b, char c) {
    return SPScore(a, b) + SPScore(b, c) + SPScore(a, c);
}

void Advance(
    const int* prev, int* cur,
    const char* u, const char* v,
    int i,
    int y0, int y1,
    int step,
    int yb0, int yb1) {
    for (int j = y0; j != y1; j += step) {
        bool assigned = false;
        int opt = 0, s, y;

        y = j - step;
        if (y >= yb0 && y < yb1) {
            s = cur[y] + SPScore('_', v[j]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        if (prev) {
            y = j;
            if (y >= yb0 && y < yb1) {
                s = prev[y] + SPScore(u[i], '_');
                if (!assigned || s > opt) {
                    assigned = true;
                    opt = s;
                }
            }

            y = j - step;
            if (y >= yb0 && y < yb1) {
                s = prev[y] + SPScore(u[i], v[j]);
                if (!assigned || s > opt) {
                    assigned = true;
                    opt = s;
                }
            }
        } // if (prev)

        if (assigned) {
            cur[j] = opt;
        }
    }
}

void Initialize(
    int* cur,
    const char* u, const char* v,
    int x,
    int y0, int y1, int step) {
    bool matched = false;
    for (int j = y0; j != y1; j += step) {
        if (matched || u[x] == v[j]) {
            cur[j] = SCORE_MATCH + SCORE_INDEL * abs(j - y0);
            matched = true;
        }
        else
            cur[j] = SCORE_MISMATCH + SCORE_INDEL * abs(j - y0);
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
    if (x1 - x0 <= 2)
        return;

    int m = (x0 + x1) / 2;

    Initialize(ctx->left, ctx->u, ctx->v, x0, y0, y1, 1);
    for (int i = 1; i <= m; ++i) {
        Advance(
            ctx->left, ctx->prev, ctx->u, ctx->v,
            i,
            y0, y1,
            1,
            y0, y1);
        swap(ctx->left, ctx->prev);
    }

    Initialize(ctx->right, ctx->u, ctx->v, x1 - 1, y1 - 1, y0 - 1, -1);
    for (int i = x1 - 2; i > m; --i) {
        Advance(
            ctx->right, ctx->prev, ctx->u, ctx->v,
            i,
            y1 - 1, y0 - 1,
            -1,
            y0, y1);
        swap(ctx->right, ctx->prev);
    }

    int k = y0, l = ctx->left[y0] + ctx->right[y0];
    for (int j = y0 + 1; j < y1; ++j) {
        int t;
        if (j == y1 - 1)
            t = (x1 - m - 1) * SCORE_INDEL;
        else
            t = ctx->right[j + 1];

        int n = ctx->left[j] + t;
        if (n > l) {
            l = n;
            k = j;
        }
    }

    ctx->path[m] = k;

    Solve(ctx, x0, y0, m + 1, k + 1);
    Solve(ctx, m, k, x1, y1);
}

int main(int argc, char* argv[]) {
    char str0[MAXSTRLEN], str1[MAXSTRLEN];

    scanf("%s%s", str0, str1);
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

    // Construct the first and the last alignment.
    ctx.path[0] = 0;
    for (int j = 0; j < ctx.path[1]; ++j) {
        if (v[j] == u[0]) {
            ctx.path[0] = j;
            break;
        }
    }

    ctx.path[lu - 1] = lv - 1;
    for (int j = ctx.path[lu - 1] + 1; j < lv; ++j) {
        if (v[j] == u[lu - 1]) {
            ctx.path[lu - 1] = j;
            break;
        }
    }

    int total = 0, length = 0, matches = 0;
    for (int i = 0, b = 0; i < lu; ++i) {
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
    printf("score: %d\nlength: %d\nmatches: %d\n", total, length, matches);

    for (int i = 0, b = 0; i < lu; ++i) {
        while (b < ctx.path[i]) {
            printf("_ ");
            ++b;
        }
        b = ctx.path[i] + 1;

        printf("%c ", u[i]);
    }
    printf("\n");

    for (int i = 0, b = 0; i < lu; ++i) {
        while (b <= ctx.path[i]) {
            printf("%c ", v[b]);
            ++b;
        }
    }
    printf("\n");


    return 0;
}
