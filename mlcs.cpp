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

    int score[2 * MAXSTRLEN];
    int* prev = score;
    int* cur = score + MAXSTRLEN;
    
    // int step = -1;
    // int y0 = lv - 1, y1 = -1;
    int yb0 = 0, yb1 = lv;
    int step = 1;
    int y0 = 0, y1 = lv;

    // Initialize(prev, u, v, lu - 1, y0, y1, step);
    Initialize(prev, u, v, 0, y0, y1, step);
    // for (int i = lu - 2; i >= 0; --i) {
    for (int i = 1; i < lu; ++i) {
        Advance(
            prev, cur, u, v,
            i,
            y0, y1,
            step,
            yb0, yb1);
        swap(prev, cur);
    }

    // printf("%d\n", prev[0]);
    printf("%d\n", prev[lv - 1]);

    return 0;
}