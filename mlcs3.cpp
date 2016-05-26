#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>

using std::max;
using std::abs;
using std::swap;
using std::make_pair;
using std::pair;
using std::vector;
using std::string;

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

    int* link0u;
    int* link0vw;
    int* link1u;
    int* link1vw;

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

// pair<score, parent>
pair<int, pair<int, int>> Score(
    const Context* ctx,
    const int* prev, const int* cur,
    int i, int j, int k,
    int step,
    int yb0, int yb1,
    int zb0, int zb1) {
    int dimy = ctx->dimv, dimz = ctx->dimw;
    bool assigned = false;
    int opt = 0, y, z;
    int px, pyz;

    y = j - step;
    z = k;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', ctx->v[j], '_');
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;

            px = i;
            pyz = y * dimz + z;
        }
    }

    y = j;
    z = k - step;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', '_', ctx->w[k]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;

            px = i;
            pyz = y * dimz + z;
        }
    }

    y = j - step;
    z = k - step;
    if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
        int s = cur[y * dimz + z] + SPScore('_', ctx->v[j], ctx->w[k]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;

            px = i;
            pyz = y * dimz + z;
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

                px = i - 1;
                pyz = y * dimz + z;
            }
        }

        y = j - step;
        z = k;
        if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
            int s = prev[y * dimz + z] + SPScore(ctx->u[i], ctx->v[j], '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;

                px = i - 1;
                pyz = y * dimz + z;
            }
        }

        y = j;
        z = k - step;
        if (y >= yb0 && y < yb1 && z >= zb0 && z < zb1) {
            int s = prev[y * dimz + z] + SPScore(ctx->u[i], '_', ctx->w[k]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;

                px = i - 1;
                pyz = y * dimz + z;
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

                px = i - 1;
                pyz = y * dimz + z;
            }
        }
    } // if (prev)

    return make_pair(opt, make_pair(px, pyz));
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
            auto r = Score(ctx, prev, cur, i, j, k, step, yb0, yb1, zb0, zb1);
            cur[j * dimz + k] = r.first;
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
            // u[0] is '_'.
            auto r =
                Score(ctx, nullptr, cur, 0, j, k, step, yb0, yb1, zb0, zb1);
            cur[j * dimz + k] = r.first;
        }
    }
}

void Step(
    const Context* ctx,
    int* linku, int* linkvw,
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
            auto r = Score(ctx, prev, cur, i, j, k, step, yb0, yb1, zb0, zb1);
            int p = j * dimz + k;
            cur[p] = r.first;

            linku[p] = r.second.first;
            linkvw[p] = r.second.second;
        }
    }
}

void Initialize(
    const Context* ctx,
    int* linku, int* linkvw,
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

        linku[y0 * dimz + k] = 0;
        linkvw[y0 * dimz + k] = y0 * dimz + k - 1;
    }

    for (int j = y0 + 1; j != y1; j += step) {
        for (int k = z0; k != z1; k += step) {
            // u[0] is '_'.
            auto r =
                Score(ctx, nullptr, cur, 0, j, k, step, yb0, yb1, zb0, zb1);
            int p = j * dimz + k;
            cur[p] = r.first;

            linku[p] = r.second.first;
            linkvw[p] = r.second.second;
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

    int buffer[(MAXSTRLEN * 7 + 2) * MAXSTRLEN];

    Context ctx;
    ctx.prev = buffer;
    ctx.left = ctx.prev + MAXSTRLEN * MAXSTRLEN;
    ctx.right = ctx.left + MAXSTRLEN * MAXSTRLEN;
    ctx.link0u = ctx.right + MAXSTRLEN * MAXSTRLEN;
    ctx.link0vw = ctx.link0u + MAXSTRLEN * MAXSTRLEN;
    ctx.link1u = ctx.link0vw + MAXSTRLEN * MAXSTRLEN;
    ctx.link1vw = ctx.link1u + MAXSTRLEN * MAXSTRLEN;

    ctx.pathv = ctx.link1vw + MAXSTRLEN * MAXSTRLEN;
    ctx.pathw = ctx.pathv + MAXSTRLEN;

    ctx.u = u;
    ctx.v = v;
    ctx.w = w;

    ctx.dimv = lv;
    ctx.dimw = lw;

    //--------------------------------------------------------------------------

    Solve(&ctx, 0, 0, 0, lu, lv, lw);

#if 0
    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.pathv[i]);
    }
    printf("\n");
    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.pathw[i]);
    }
    printf("\n");
#endif
    
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
    vector<pair<int, int>> path, seg;
    pair<int, int> last = make_pair(0, 0);
    path.push_back(last);
    Initialize(
        &ctx,
        ctx.link0u, ctx.link0vw,
        ctx.left,
        0, 0, lv, lw,
        1,
        0, lv,
        0, lw);
    for (int i = 1; i < lu; ++i) {
        Step(&ctx,
            ctx.link1u, ctx.link1vw,
            ctx.left, ctx.prev,
            i,
            0, 0,
            lv, lw,
            1,
            0, lv,
            0, lw);

        seg.clear();
        int* linku = ctx.link1u;
        int* linkvw = ctx.link1vw;
        pair<int, int> c = make_pair(i, ctx.pathv[i] * lw + ctx.pathw[i]);
        while (c != last) {
            seg.push_back(c);
            int px = linku[c.second];
            int py = linkvw[c.second];
            if (px != i) {
                linku = ctx.link0u;
                linkvw = ctx.link0vw;
            }
            c = make_pair(px, py);
        }

        for (int t = seg.size() - 1; t >= 0; --t) {
            path.push_back(seg[t]);
        }
        last = make_pair(i, ctx.pathv[i] * lw + ctx.pathw[i]);

        swap(ctx.left, ctx.prev);
        swap(ctx.link0u, ctx.link1u);
        swap(ctx.link0vw, ctx.link1vw);
    }
    seg.clear();
    pair<int, int> c = make_pair(lu - 1, (lv - 1) * lw + lw - 1);
    while (c != last) {
        seg.push_back(c);
        int px = ctx.link0u[c.second];
        int py = ctx.link0vw[c.second];
        c = make_pair(px, py);
    }
    for (int t = seg.size() - 1; t >= 0; --t) {
        path.push_back(seg[t]);
    }

#if 0
    for (int i = 0; i < path.size(); ++i) {
        printf("%d %d %d\n", path[i].first,
            path[i].second / lw, path[i].second % lw);
    }
#endif

    int total = 0, length = path.size() - 1, matches = 0;
    string au, av, aw;
    for (int q = 1, i = 1, j = 1, k = 1; q < path.size(); ++q) {
        char a = '_';
        if (path[q].first != path[q - 1].first) {
            a = u[i];
            ++i;
        }

        int v1 = path[q].second / lw;
        int v0 = path[q - 1].second / lw;
        char b = '_';
        if (v1 != v0) {
            b = v[j];
            ++j;
        }

        int w1 = path[q].second % lw;
        int w0 = path[q - 1].second % lw;
        char c = '_';
        if (w1 != w0) {
            c = w[k];
            ++k;
        }

        total += SPScore(a, b, c);
        if (a == b && b == c && a != '_') {
            ++matches;
        }

        au.push_back(a); au.push_back(' ');
        av.push_back(b); av.push_back(' ');
        aw.push_back(c); aw.push_back(' ');
    }

    printf("score: %d\nlength: %d\nmatches: %d\n", total, length, matches);

    printf("%s\n%s\n%s\n", au.c_str(), av.c_str(), aw.c_str());

    return 0;
}