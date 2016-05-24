#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstring>

using std::max;
using std::swap;

const int MAXSTRLEN = 128;

int SPScore(char a, char b) {
    if (a == b && a != '_') {
        return 5;
    }
    else if (a != b && a != '_' && b != '_') {
        return -4;
    }
    else if (a != '_' || b != '_') {
        return -8;
    }
    else
        return 0;
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

int main(int argc, char* argv[]) {
    char u[MAXSTRLEN], v[MAXSTRLEN];

    scanf("%s%s", u, v);
    int lu = strlen(u), lv = strlen(v);

    int score[2 * MAXSTRLEN];
    int* prev = score;
    int* cur = score + MAXSTRLEN;

    cur[0] = max(SPScore('_', '_'),
        max(SPScore('_', v[0]),
            max(SPScore(u[0], '_'), SPScore(u[0], v[0]))));
    
    int step = 1;
    int y0 = 0, y1 = lv;
    int yb0 = 0, yb1 = lv;

    Advance(
        nullptr, cur, u, v,
        0,
        y0, y1,
        step,
        yb0, yb1);
    swap(prev, cur);
    for (int i = 1; i < lu; ++i) {
        Advance(
            prev, cur, u, v,
            i,
            y0, y1,
            step,
            yb0, yb1);
        swap(prev, cur);
    }

    printf("%d\n", prev[lv - 1]);    

    return 0;
}