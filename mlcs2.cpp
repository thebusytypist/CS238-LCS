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
        cur[j] = SCORE_INDEL * abs(j - y0);
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

// Return the maximum score.
int Solve(Context* ctx, int x0, int y0, int x1, int y1) {
    if (x1 - x0 <= 1)
        return 0;

    // Handle the degenerate case.
    if (y0 == y1) {
        for (int i = x0 + 1; i < x1; ++i) {
            ctx->path[i] = 0;
        }
        // We do not care the return value for the non-top case.
        return 0;
    }

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

    return l;
}

void Trace(vector<pair<int, int>>& path,
    const int* prev, const int* cur,
    const char* u, const char* v,
    int x0, int y0,
    int x1, int y1) {
    // pair<parent, pair<x, y>>
    typedef pair<int, pair<int, int>> Node;
    vector<Node> Q;
    Q.push_back(make_pair(-1, make_pair(x1, y1)));
    int front = 0;

    while (front < Q.size()) {
        auto h = Q[front];
        int x = h.second.first;
        int y = h.second.second;

        if (x == x0 && y == y0)
            break;

        int best = 0;
        bool assigned = false;

        bool s0valid = x > x0;
        int s0 = 0;
        if (s0valid) {
            s0 = prev[y] + SPScore(u[x], '_');
            best = s0 > best || !assigned ? s0 : best;
            assigned = true;
        }

        bool s1valid = x > x0 && y > y0;
        int s1 = 0;
        if (s1valid) {
            s1 = prev[y - 1] + SPScore(u[x], v[y]);
            best = s1 > best || !assigned ? s1 : best;
            assigned = true;
        }

        bool s2valid = y > y0;
        int s2 = 0;
        if (s2valid) {
            const int* p = x == x1 ? cur : prev;
            s2 = p[y - 1] + SPScore('_', v[y]);
            best = s2 > best || !assigned ? s2 : best;
            assigned = true;
        }

        if (s0valid && s0 == best) {
            Q.push_back(make_pair(front, make_pair(x - 1, y)));
        }
        if (s1valid && s1 == best) {
            Q.push_back(make_pair(front, make_pair(x - 1, y - 1)));
        }
        if (s2valid && s2 == best) {
            Q.push_back(make_pair(front, make_pair(x, y - 1)));
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

    int score = Solve(&ctx, 0, 0, lu, lv);

    printf("score: %d\n", score);

#if 0
    for (int i = 1; i < lu; ++i) {
        printf("%d ", ctx.path[i]);
    }
    printf("\n");
#endif

    //--------------------------------------------------------------------------
    vector<pair<int, int>> path;
    path.push_back(make_pair(0, 0));

    Initialize(ctx.left, ctx.u, ctx.v, 0, lv, 1);
    for (int i = 1; i < lu; ++i) {
        Step(
            ctx.left, ctx.prev, ctx.u, ctx.v,
            i,
            0, lv,
            1,
            0, lv);
        swap(ctx.left, ctx.prev);
        int last = path[path.size() - 1].second;
        Trace(path, ctx.prev, ctx.left, u, v,
            i - 1, last,
            i, ctx.path[i]);
    }
    int last = path[path.size() - 1].second;
    if (last != lv - 1) {
        Trace(path, nullptr, ctx.left, u, v,
            lu - 1, last,
            lu - 1, lv - 1);
    }
    printf("reference score: %d\n", ctx.left[lv - 1]);

    //--------------------------------------------------------------------------
#if 0
    for (int t = 0; t < path.size(); ++t) {
        printf("(%d %d)\n", path[t].first, path[t].second);
    }
#endif

    int length = path.size() - 1, matches = 0;
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
        if (a == b && a != '_') {
            ++matches;
        }

        au.push_back(a); au.push_back(' ');
        av.push_back(b); av.push_back(' ');
    }

    printf("length: %d\nmatches: %d\n", length, matches);

    printf("%s\n%s\n", au.c_str(), av.c_str());

    return 0;
}
