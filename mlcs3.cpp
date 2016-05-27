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

#ifndef TEST
const int MAXSTRLEN = 128;
#endif

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
    int dimz = ctx->dimw;
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

int Score2(
    const int* prev, const int* cur,
    const char* u, const char* v,
    int i, int j,
    int step,
    int yb0, int yb1) {
    // The last parameter of SPScore is '_' by default.
    bool assigned = false;
    int opt = 0, y;

    y = j - step;
    if (y >= yb0 && y < yb1) {
        int s = cur[y] + SPScore('_', v[j], '_');
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
        }
    }

    if (prev) {
        y = j;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], '_', '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
            }
        }

        y = j - step;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], v[j], '_');
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
    for (int k = z0 + step; k != z1; k += step) {
        cur[y0 * dimz + k] = 2 * SCORE_INDEL * abs(k - z0);
    }

    for (int j = y0 + step; j != y1; j += step) {
        for (int k = z0; k != z1; k += step) {
            // u[0] is '_'.
            cur[j * dimz + k] =
                Score(ctx, nullptr, cur, 0, j, k, step, yb0, yb1, zb0, zb1);
        }
    }
}

void Step2(
    const int* prev, int* cur,
    const char* u, const char* v,
    int i,
    int y0, int y1,
    int step,
    int yb0, int yb1) {
    for (int j = y0; j != y1; j += step) {
        cur[j] = Score2(prev, cur, u, v, i, j, step, yb0, yb1);
    }
}

void Initialize2(
    int* cur,
    const char* u, const char* v,
    int y0, int y1, int step) {
    cur[y0] = SCORE_SPACES;
    for (int j = y0 + step; j != y1; j += step) {
        // The score is for a match of _, _, *.
        cur[j] = 2 * SCORE_INDEL * abs(j - y0);
    }
}

int Solve2(
    int* path,
    int* prev, int* left, int* right,
    const char* u, const char* v,
    int x0, int y0, int x1, int y1) {
    if (x1 - x0 <= 1)
        return 0;

    // Handle the degenerate case.
    if (y0 == y1) {
        for (int i = x0 + 1; i < x1; ++i) {
            path[i] = y0;
        }
        // We do not care the return value for the non-top case.
        return 0;
    }

    int m = (x0 + x1) / 2;

    Initialize2(left, u, v, y0, y1, 1);
    for (int i = x0 + 1; i <= m; ++i) {
        Step2(
            left, prev, u, v,
            i,
            y0, y1,
            1,
            y0, y1);
        swap(left, prev);
    }

    Initialize2(right, u, v, y1, y0, -1);
    for (int i = x1 - 1; i > m; --i) {
        Step2(
            right, prev, u, v,
            i,
            y1, y0,
            -1,
            y0 + 1, y1 + 1);
        swap(right, prev);
    }

    int k = y0, l = left[y0] + right[y0 + 1];
    for (int j = y0 + 1; j < y1; ++j) {
        int n = left[j] + right[j + 1];
        if (n > l) {
            l = n;
            k = j;
        }
    }

    path[m] = k;

    Solve2(path, prev, left, right, u, v, x0, y0, m, k);
    Solve2(path, prev, left, right, u, v, m, k, x1, y1);

    return l;
}

int Solve(Context* ctx, int x0, int y0, int z0, int x1, int y1, int z1) {
    if (x1 - x0 <= 1)
        return 0;

    // Handle the degenerate cases.
    if (y0 == y1) {
        for (int i = x0 + 1; i < x1; ++i) {
            ctx->pathv[i] = y0;
        }

        Solve2(ctx->pathw,
            ctx->prev, ctx->left, ctx->right,
            ctx->u, ctx->w,
            x0, z0, x1, z1);
        return 0;
    }
    else if (z0 == z1) {
        for (int i = x0 + 1; i < x1; ++i) {
            ctx->pathw[i] = z0;
        }

        Solve2(ctx->pathv,
            ctx->prev, ctx->left, ctx->right,
            ctx->u, ctx->v,
            x0, y0, x1, y1);
        return 0;
    }

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

    return l;
}

typedef pair<int, pair<int, int>> Position;

void Trace(
    const Context* ctx,
    vector<Position>& path,
    const int* prev, const int* cur,
    int x0, int y0, int z0,
    int x1, int y1, int z1) {
    int dimz = ctx->dimw;
    const char* u = ctx->u;
    const char* v = ctx->v;
    const char* w = ctx->w;
    // pair<parent, pair<x, y, z>>
    typedef pair<int, Position> Node;
    vector<Node> Q;
    Q.push_back(make_pair(-1, make_pair(x1, make_pair(y1, z1))));
    int front = 0;

    while (front < Q.size()) {
        auto h = Q[front];
        int x = h.second.first;
        int y = h.second.second.first;
        int z = h.second.second.second;

        if (x == x0 && y == y0 && z == z0)
            break;

        int best = 0;
        bool assigned = false;
        const int* p = x == x1 ? cur : prev;

        bool s0valid = y > y0;
        int s0 = 0;
        if (s0valid) {
            s0 = p[(y - 1) * dimz + z] + SPScore('_', v[y], '_');
            best = s0 > best || !assigned ? s0 : best;
            assigned = true;
        }

        bool s1valid = y > y0 && z > z0;
        int s1 = 0;
        if (s1valid) {
            s1 = p[(y - 1) * dimz + z - 1] + SPScore('_', v[y], w[z]);
            best = s1 > best || !assigned ? s1 : best;
            assigned = true;
        }

        bool s2valid = z > z0;
        int s2 = 0;
        if (s2valid) {
            s2 = p[y * dimz + z - 1] + SPScore('_', '_', w[z]);
            best = s2 > best || !assigned ? s2 : best;
            assigned = true;
        }

        //----------------------------------------------------------------------

        bool prevvalid = x > x0;

        bool s3valid = prevvalid;
        int s3 = 0;
        if (s3valid) {
            s3 = prev[y * dimz + z] + SPScore(u[x], '_', '_');
            best = s3 > best || !assigned ? s3 : best;
            assigned = true;
        }

        bool s4valid = prevvalid && y > y0;
        int s4 = 0;
        if (s4valid) {
            s4 = prev[(y - 1) * dimz + z] + SPScore(u[x], v[y], '_');
            best = s4 > best || !assigned ? s4 : best;
            assigned = true;
        }

        bool s5valid = prevvalid && y > y0 && z > z0;
        int s5 = 0;
        if (s5valid) {
            s5 = prev[(y - 1) * dimz + z - 1] + SPScore(u[x], v[y], w[z]);
            best = s5 > best || !assigned ? s5 : best;
            assigned = true;
        }

        bool s6valid = prevvalid && z > z0;
        int s6 = 0;
        if (s6valid) {
            s6 = prev[y * dimz + z - 1] + SPScore(u[x], '_', w[z]);
            best = s6 > best || !assigned ? s6 : best;
            assigned = true;
        }

        if (s0valid && s0 == best) {
            Q.push_back(make_pair(front, make_pair(x, make_pair(y - 1, z))));
        }
        if (s1valid && s1 == best) {
            Q.push_back(make_pair(front, make_pair(x,
                make_pair(y - 1, z - 1))));
        }
        if (s2valid && s2 == best) {
            Q.push_back(make_pair(front, make_pair(x, make_pair(y, z - 1))));
        }

        if (s3valid && s3 == best) {
            Q.push_back(make_pair(front,
                make_pair(x - 1, make_pair(y, z))));
        }
        if (s4valid && s4 == best) {
            Q.push_back(make_pair(front,
                make_pair(x - 1, make_pair(y - 1, z))));
        }
        if (s5valid && s5 == best) {
            Q.push_back(make_pair(front,
                make_pair(x - 1, make_pair(y - 1, z - 1))));
        }
        if (s6valid && s6 == best) {
            Q.push_back(make_pair(front,
                make_pair(x - 1, make_pair(y, z - 1))));
        }

        // Pop the current node.
        ++front;
    }

    front = Q[front].first;
    while (front != -1) {
        path.push_back(Q[front].second);
        front = Q[front].first;
    }
}

int main() {
#ifndef TEST
    char str0[MAXSTRLEN], str1[MAXSTRLEN], str2[MAXSTRLEN];
    scanf("%s%s%s", str0 + 1, str1 + 1, str2 + 1);
#else
    int sz0, sz1, sz2;
    scanf("%d", &sz0);
    char* str0 = (char*)malloc(sizeof(char) * (sz0 + 2));
    scanf("%s", str0 + 1);

    scanf("%d", &sz1);
    char* str1 = (char*)malloc(sizeof(char) * (sz1 + 2));
    scanf("%s", str1 + 1);

    scanf("%d", &sz2);
    char* str2 = (char*)malloc(sizeof(char) * (sz2 + 2));
    scanf("%s", str2 + 1);
#endif

    str0[0] = str1[0] = str2[0] = '_';

    const char* u = str0;
    const char* v = str1;
    const char* w = str2;

#ifndef TEST
    int lu = strlen(u), lv = strlen(v), lw = strlen(w);
#else
    int lu = sz0 + 1, lv = sz1 + 1, lw = sz2 + 1;
#endif

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

    // printf("%s %s %s\n", u + 1, v + 1, w + 1);

#ifndef TEST
    int buffer[(MAXSTRLEN * 3 + 2) * MAXSTRLEN];
#else
    int rowUnit = (lu + 1);
    int pageUnit = (lv + 1) * (lw + 1);
    int* buffer = (int*)malloc((3 * pageUnit + 2 * rowUnit) * sizeof(int));
#endif

    Context ctx;
    ctx.prev = buffer;
#ifndef TEST
    ctx.left = ctx.prev + MAXSTRLEN * MAXSTRLEN;
    ctx.right = ctx.left + MAXSTRLEN * MAXSTRLEN;

    ctx.pathv = ctx.right + MAXSTRLEN * MAXSTRLEN;
    ctx.pathw = ctx.pathv + MAXSTRLEN;
#else
    ctx.left = ctx.prev + pageUnit;
    ctx.right = ctx.left + pageUnit;

    ctx.pathv = ctx.right + pageUnit;
    ctx.pathw = ctx.pathv + rowUnit;
#endif

    ctx.u = u;
    ctx.v = v;
    ctx.w = w;

    ctx.dimv = lv;
    ctx.dimw = lw;

    //--------------------------------------------------------------------------

    int score = Solve(&ctx, 0, 0, 0, lu, lv, lw);

    printf("score: %d\n", score);

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
    vector<Position> path;
    path.push_back(make_pair(0, make_pair(0, 0)));

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

        auto last = path[path.size() - 1].second;
        Trace(&ctx, path,
            ctx.prev, ctx.left,
            i - 1, last.first, last.second,
            i, ctx.pathv[i], ctx.pathw[i]);
    }
    auto last = path[path.size() - 1].second;
    if (last.first != lv - 1 || last.second != lw - 1) {
        Trace(&ctx, path,
            nullptr, ctx.left,
            lu - 1, last.first, last.second,
            lu - 1, lv - 1, lw - 1);
    }

    printf("reference score: %d\n", ctx.left[(lv - 1) * lw + lw - 1]);

    //--------------------------------------------------------------------------
#if 0
    for (auto& q : path) {
        printf("(%d %d %d)\n", q.first, q.second.first, q.second.first);
    }
#endif

    int length = path.size() - 1, matches = 0;
    string au, av, aw;
    for (int q = 1, i = 1, j = 1, k = 1; q < path.size(); ++q) {
        char a = '_';
        if (path[q].first != path[q - 1].first) {
            a = u[i];
            ++i;
        }

        int v1 = path[q].second.first;
        int v0 = path[q - 1].second.first;
        char b = '_';
        if (v1 != v0) {
            b = v[j];
            ++j;
        }

        int w1 = path[q].second.second;
        int w0 = path[q - 1].second.second;
        char c = '_';
        if (w1 != w0) {
            c = w[k];
            ++k;
        }

        if (a == b && b == c && a != '_') {
            ++matches;
        }

#if 0
        au.push_back(a); au.push_back(' ');
        av.push_back(b); av.push_back(' ');
        aw.push_back(c); aw.push_back(' ');
#endif
    }

    printf("length: %d\nmatches: %d\n", length, matches);

#if 0
    printf("%s\n%s\n%s\n", au.c_str(), av.c_str(), aw.c_str());
#endif
    return 0;
}