#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>

using std::min;
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

// pair<score, parent>
pair<int, pair<int, int>> Score(
    const int* prev, const int* cur,
    const char* u, const char* v,
    int i, int j,
    int step,
    int yb0, int yb1) {
    bool assigned = false;
    int opt = 0, y, py, px;

    y = j - step;
    if (y >= yb0 && y < yb1) {
        int s = cur[y] + SPScore('_', v[j]);
        if (!assigned || s > opt) {
            assigned = true;
            opt = s;
            py = y;
            px = i;
        }
    }

    if (prev) {
        y = j;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], '_');
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
                py = y;
                px = i - 1;
            }
        }

        y = j - step;
        if (y >= yb0 && y < yb1) {
            int s = prev[y] + SPScore(u[i], v[j]);
            if (!assigned || s > opt) {
                assigned = true;
                opt = s;
                py = y;
                px = i - 1;
            }
        }
    } // if (prev)

    return make_pair(opt, make_pair(px, py));
}

void Step(
    const int* prev, int* cur,
    const char* u, const char* v,
    int i,
    int y0, int y1,
    int step,
    int yb0, int yb1) {
    for (int j = y0; j != y1; j += step) {
        auto r = Score(prev, cur, u, v, i, j, step, yb0, yb1);
        cur[j] = r.first;
    }
}

void Initialize(
    int* cur,
    const char* u, const char* v,
    int y0, int y1, int step) {
    cur[y0] = SCORE_SPACES;
    for (int j = y0 + 1; j != y1; j += step) {
        cur[j] = SCORE_INDEL * abs(j - y0);
    }
}

void Step(
    int* linku, int* linkv,
    const int* prev, int* cur,
    const char* u, const char* v,
    int i,
    int y0, int y1,
    int step,
    int yb0, int yb1) {
    for (int j = y0; j != y1; j += step) {
        auto r = Score(prev, cur, u, v, i, j, step, yb0, yb1);
        cur[j] = r.first;

        linku[j] = r.second.first;
        linkv[j] = r.second.second;
    }
}

void Initialize(
    int* linku, int* linkv,
    int* cur,
    const char* u, const char* v,
    int y0, int y1, int step) {
    cur[y0] = SCORE_SPACES;
    for (int j = y0 + 1; j != y1; j += step) {
        cur[j] = SCORE_INDEL * abs(j - y0);

        linku[j] = 0;
        linkv[j] = j - 1;
    }
}

struct Context {
    int* prev;
    int* left;
    int* right;
    int* path;

    int* link0u;
    int* link0v;
    int* link1u;
    int* link1v;

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

    Initialize(ctx->right, ctx->u, ctx->v, y1, y0, -1);
    for (int i = x1 - 1; i > m; --i) {
        Step(
            ctx->right, ctx->prev, ctx->u, ctx->v,
            i,
            y1, y0,
            -1,
            y0 + 1, y1 + 1);
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

    int buffer[8 * MAXSTRLEN];

    Context ctx;

    ctx.prev = buffer;
    ctx.left = ctx.prev + MAXSTRLEN;
    ctx.right = ctx.left + MAXSTRLEN;
    ctx.path = ctx.right + MAXSTRLEN;
    ctx.link0u = ctx.path + MAXSTRLEN;
    ctx.link1u = ctx.link0u + MAXSTRLEN;
    ctx.link0v = ctx.link1u + MAXSTRLEN;
    ctx.link1v = ctx.link0v + MAXSTRLEN;
    ctx.u = u;
    ctx.v = v;

    Solve(&ctx, 0, 0, lu, lv);

#if 0
    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.path[i]);
    }
    printf("\n");
#endif

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
    ctx.path[0] = 0;
    vector<pair<int, int>> path, seg;
    pair<int, int> last = make_pair(0, 0);
    path.push_back(last);
    Initialize(ctx.link0u, ctx.link0v, ctx.left, ctx.u, ctx.v, 0, lv, 1);
    for (int i = 1; i < lu; ++i) {
        Step(
            ctx.link1u, ctx.link1v,
            ctx.left, ctx.prev, ctx.u, ctx.v,
            i,
            0, lv,
            1,
            0, lv);

        seg.clear();
        int* linku = ctx.link1u;
        int* linkv = ctx.link1v;
        pair<int, int> c = make_pair(i, ctx.path[i]);
        while (c != last) {
            seg.push_back(c);
            int px = linku[c.second];
            int py = linkv[c.second];
            if (px != i) {
                linku = ctx.link0u;
                linkv = ctx.link0v;
            }
            c = make_pair(px, py);
        }

        for (int t = seg.size() - 1; t >= 0; --t) {
            path.push_back(seg[t]);
        }
        last = make_pair(i, ctx.path[i]);

        swap(ctx.left, ctx.prev);
        swap(ctx.link0u, ctx.link1u);
        swap(ctx.link0v, ctx.link1v);
    }
    seg.clear();
    pair<int, int> c = make_pair(lu - 1, lv - 1);
    while (c != last) {
        seg.push_back(c);
        int px = ctx.link0u[c.second];
        int py = ctx.link0v[c.second];
        c = make_pair(px, py);
    }
    for (int t = seg.size() - 1; t >= 0; --t) {
        path.push_back(seg[t]);
    }

#if 0
    for (int i = 0; i < path.size(); ++i) {
        printf("%d %d\n", path[i].first, path[i].second);
    }
#endif

    int total = 0, length = path.size() - 1, matches = 0;
    string au, av;
    for (int q = 1, i = 1, j = 1; q < path.size(); ++q) {
        char a = '_';
        if (path[q].first != path[q - 1].first) {
            a = u[i];
            ++i;
        }
        char b = '_';
        if (path[q].second != path[q - 1].second) {
            b = v[j];
            ++j;
        }
        total += SPScore(a, b);
        if (a == b && a != '_') {
            ++matches;
        }

        au.push_back(a); au.push_back(' ');
        av.push_back(b); av.push_back(' ');
    }

    printf("score: %d\nlength: %d\nmatches: %d\n", total, length, matches);

    printf("%s\n%s\n", au.c_str(), av.c_str());

    return 0;
}
