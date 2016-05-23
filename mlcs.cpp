#include <algorithm>
#include <cstdio>
#include <cstring>
using std::max;

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

int main(int argc, char* argv[]) {
    char u[MAXSTRLEN], v[MAXSTRLEN];

    scanf("%s%s", u, v);
    int lu = strlen(u), lv = strlen(v);

    int score[MAXSTRLEN][MAXSTRLEN];
    score[0][0] = max(SPScore('_', '_'),
        max(SPScore('_', v[0]),
            max(SPScore(u[0], '_'), SPScore(u[0], v[0]))));
    
    for (int i = 0; i < lu; ++i) {
        for (int j = 0; j < lv; ++j) {
            if (i == 0 && j == 0)
                continue;

            bool assigned = false;
            int opt = 0;

            int x = i - 1, y = j - 1;
            if (x >= 0 && y >= 0) {
                int s = score[x][y] + SPScore(u[i], v[j]);
                if (!assigned || s > opt) {
                    opt = s;
                    assigned = true;
                }
            }

            x = i - 1; y = j;
            if (x >= 0 && y >= 0) {
                int s = score[x][y] + SPScore(u[i], '_');
                if (!assigned || s > opt) {
                    opt = s;
                    assigned = true;
                }
            }

            x = i; y = j - 1;
            if (x >= 0 && y >= 0) {
                int s = score[x][y] + SPScore('_', v[j]);
                if (!assigned || s > opt) {
                    opt = s;
                    assigned = true;
                }
            }

            score[i][j] = opt;
        }
    }

    for (int j = 0; j < lv; ++j) {
        for (int i = 0; i < lu; ++i) {
            printf("%d ", score[i][j]);
        }
        printf("\n");
    }

    return 0;
}