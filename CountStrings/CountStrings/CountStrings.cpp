// CountStrings.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>

using namespace std;

#define MOD 1000000007

struct node;
struct edge;
struct graph;
struct exp_state;
typedef struct node Node;
typedef struct edge Edge;
typedef struct graph Graph;
typedef exp_state ExpState;

template <typename T>
class MyVec {
    T* buf;
    int capa;
    int idx;
public:
    MyVec(int capa = 1) : capa(capa), idx(0) {
        buf = new T[capa];
    }
    ~MyVec() {
        delete[] buf;
    }
    void push_back(const T& c) {
        if (idx == capa) {
            T* tmp = new T[capa << 1];
            for (int i = 0; i < capa; ++i) tmp[i] = buf[i];
            delete[] buf;
            buf = tmp;
            capa <<= 1;
        }
        buf[idx++] = c;
    }
    int size() const { return idx; }
    T& operator[](int i) { return buf[i]; }
    void clear() { idx = 0; }
};

class IntSet {
    int len;
    int* data;
    int set_size;
public:
    // *********************** iterator ************************************
    class iterator {
        int* data;
        int pos;
        int capa;
    public:
        iterator(int* data=nullptr, int pos=-1, int capa=0) : data(data), pos(pos), capa(capa) {}
        bool operator==(const iterator& other) const {
            return data == other.data && pos == other.pos && capa == other.capa;
        }
        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }
        iterator& operator++() {
            while (data[++pos] == -1 && pos < capa);
            return *this;
        }
        int& operator*() {
            return data[pos];
        }
    };
    // *********************** iterator ************************************
    IntSet(int sz=4) : len(sz), set_size(0) {
        data = new int[len];
        clear();
    }
    IntSet(const IntSet& other) : len(other.len), set_size(other.set_size) {
        data = new int[len];
        for (int i = 0; i < len; ++i) data[i] = other.data[i];
    }
    ~IntSet() {
        delete[] data;
    }
    int size() const { return set_size; }
    void clear() {
        for (int i = 0; i < len; ++i) data[i] = -1;
        set_size = 0;
    }
    bool insert(int val) {
        if (val < len) {
            if (data[val] == -1) ++set_size;
            data[val] = val;
            return true;
        }
        else {
            int new_len = len;
            while (new_len <= val) new_len <<= 1;
            int* tmp = new int[new_len];
            for (int i = 0; i < new_len; ++i) tmp[i] = -1;
            for (int i = 0; i < len; ++i) tmp[i] = data[i];
            len = new_len;
            delete[] data;
            data = tmp;
            data[val] = val;
            ++set_size;
            return true;
        }
        return false;
    }
    void insert(iterator begin, iterator end) {
        iterator iter = begin;
        while (iter != end) {
            insert(*iter);
            ++iter;
        }
    }
    bool operator==(const IntSet& other) {
        if (len != other.len) return false;
        if (set_size != other.set_size) return false;
        for (int i = 0; i < len; ++i) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }
    IntSet& operator=(const IntSet& other) {
        if (this != &other) {
            len = other.len;
            set_size = other.set_size;
            delete[] data;
            data = new int[len];
            for (int i = 0; i < len; ++i) data[i] = other.data[i];
        }
        return *this;
    }
    iterator begin() {
        int pos = -1;
        while (data[++pos] == -1 && pos < len);
        return iterator(data, pos, len);
    }
    iterator end() {
        return iterator(data, len, len);
    }
    iterator find(int val) {
        if (val < len && data[val]!=-1) return iterator(data, val, len);
        return iterator(data, len, len);
    }
};

template <typename T>
class MyQ {
    T* buf;
    int capa;
    int head;
    int tail;
public:
    MyQ(int capa = 2) : capa(capa), head(0), tail(0) {
        buf = new T[capa];
    }
    ~MyQ() {
        delete[] buf;
    }
    void push(T& val) {
        if (size() == capa) {
            T* tmp = new T[capa << 1];
            for (int i = 0; i < capa; ++i) {
                tmp[i] = buf[(i + tail) % capa];
            }
            tail = 0;
            head = capa;
            capa <<= 1;
            delete[] buf;
            buf = tmp;
        }
        buf[head % capa] = val;
        ++head;
    }
    int size() const {
        return head - tail;
    }
    T& front() {
        return buf[tail % capa];
    }
    void pop() {
        ++tail;
    }
};

struct edge {
    char input;
    Node* to;
};

struct node {
    int idx;
    bool accept;
    MyVec<Edge> ne;
    IntSet closure;
};

struct graph {
    Node* start;
    Node* end;
};

struct exp_state {
    int n_atom;
    int n_alt;
};

#define TBL_SZ 200
Node ntbl[TBL_SZ];
int ntbl_cnt;
Node dtbl[TBL_SZ];
int dtbl_cnt;

void prt_mtx(int dim, int** a) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

int** alloc_sqr_mtx(int dim) {
    int** ans = new int* [dim];
    for (int i = 0; i < dim; ++i) {
        ans[i] = new int[dim];
        for (int j = 0; j < dim; ++j) {
            ans[i][j] = 0;
        }
    }
    return ans;
}

void free_sqr_mtx(int dim, int** mtx) {
    for (int i = 0; i < dim; ++i) {
        delete[] mtx[i];
    }
    delete[] mtx;
}


// c = a*b, where a and b are square matrices with dimension of dim
void mtx_mult_mod(int dim, int** a, int** b, int** c) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            c[i][j] = 0;
            for (int k = 0; k < dim; ++k) {
                long long tmp = (long long)a[i][k] * b[k][j];
                tmp %= MOD;
                tmp += c[i][j];
                tmp %= MOD;
                c[i][j] = (int)tmp;
            }
        }
    }
}

void mtx_cpy(int dim, int** src, int** dst) {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}


void mtx_pwr(int dim, int n, int** src, int** dst) {
    if (n == 1) {
        mtx_cpy(dim, src, dst);
    }
    else if (n == 2) {
        mtx_mult_mod(dim, src, src, dst);
    }
    else {
        int** buf = alloc_sqr_mtx(dim);
        mtx_pwr(dim, n / 2, src, buf);
        mtx_pwr(dim, 2, buf, dst);
        if (n % 2) {
            mtx_cpy(dim, dst, buf);
            mtx_mult_mod(dim, buf, src, dst);
        }
        free_sqr_mtx(dim, buf);
    }
}

string in2post(string& in) {
    string ans;
    string wrong("Illegal expression");
    ExpState stk[TBL_SZ], *stkp;
    ExpState es;
    int n_atom = 0;
    int n_alt = 0;
    int n = in.size();
    stkp = stk;
#define push(a) *stkp++ = a
#define pop() *--stkp;
    for (int i = 0; i < n; ++i) {
        switch (in[i]) {
        case '(':
            if (n_atom > 1) {
                --n_atom;
                ans.append(1, '.');
            }
            es.n_atom = n_atom;
            es.n_alt = n_alt;
            push(es);
            n_atom = n_alt = 0;
            break;
        case ')':
            if (stkp == stk) return wrong;
            if (n_atom == 0) return wrong;
            if (n_atom > 1) {
                --n_atom;
                ans.append(1, '.');
            }
            while (n_alt-- > 0) ans.append(1, '|');
            es = pop();
            n_atom = es.n_atom;
            n_alt = es.n_alt;
            ++n_atom;
            break;
        case '*':
        case '+':
        case '?':
            if (n_atom == 0) return wrong;
            ans.append(1, in[i]);
            break;
        case '|':
            if (n_atom == 0) return wrong;
            while (--n_atom > 0) ans.append(1, '.');
            ++n_alt;
            break;
        default:
            if (n_atom > 1) {
                --n_atom;
                ans.append(1, '.');
            }
            ++n_atom;
            ans.append(1, in[i]);
            break;
        }
    }
    if (stk != stkp) return wrong;
    if (n_atom > 1) {
        --n_atom;
        ans.append(1, '.');
    }
    while (n_alt-- > 0) ans.append(1, '|');
#undef push
#undef pop
    return ans;
}

Graph concat(Graph& a, Graph& b) {
    Node *br = &ntbl[ntbl_cnt++];

    Edge edge;
    edge.input = 'e';
    edge.to = b.start;
    br->ne.push_back(edge);
    edge.to = br;
    a.end->ne.push_back(edge);

    Graph ans;
    ans.start = a.start;
    ans.end = b.end;
    return ans;
}

Graph alt(Graph& a, Graph& b) {
    Node* start = &ntbl[ntbl_cnt++];
    Node* end = &ntbl[ntbl_cnt++];

    Edge edge;
    edge.input = 'e';
    edge.to = a.start;
    start->ne.push_back(edge);
    edge.to = b.start;
    start->ne.push_back(edge);
    edge.to = end;
    a.end->ne.push_back(edge);
    b.end->ne.push_back(edge);

    Graph ans;
    ans.start = start;
    ans.end = end;
    return ans;
}

Graph star(Graph& a) {
    Node* start = &ntbl[ntbl_cnt++];
    Node* end = &ntbl[ntbl_cnt++];

    Edge edge;
    edge.input = 'e';
    edge.to = a.start;
    start->ne.push_back(edge);
    edge.to = end;
    start->ne.push_back(edge);
    a.end->ne.push_back(edge);
    edge.to = a.start;
    a.end->ne.push_back(edge);

    Graph ans;
    ans.start = start;
    ans.end = end;
    return ans;
}

Graph atom(char c) {
    Node* start = &ntbl[ntbl_cnt++];
    Node* end = &ntbl[ntbl_cnt++];
    Edge edge;
    edge.input = c;
    edge.to = end;
    start->ne.push_back(edge);

    Graph ans;
    ans.start = start;
    ans.end = end;
    return ans;
}

Graph build_nfa(string& re) {
    int n = re.size();
    Graph stk[TBL_SZ], * stkp;
    stkp = stk;
#define push(x) *stkp++ = x
#define pop() *--stkp;
    Graph g1, g2;
    for (int i = 0; i < n; ++i) {
        switch (re[i]) {
        default:
            push(atom(re[i]));
            break;
        case '.':
            g1 = pop();
            g2 = pop();
            push(concat(g2, g1));
            break;
        case '|':
            g1 = pop();
            g2 = pop();
            push(alt(g2, g1));
            break;
        case '*':
            g1 = pop();
            push(star(g1));
            break;
        }
    }
    g1 = pop();
#undef push
#undef pop
    return g1;
}

void init_tbl() {
    for (int i = 0; i < TBL_SZ; ++i) {
        ntbl[i].ne.clear();
        ntbl[i].closure.clear();
        ntbl[i].accept = false;
        ntbl[i].idx = i;

        dtbl[i].ne.clear();
        dtbl[i].closure.clear();
        dtbl[i].accept = false;
        dtbl[i].idx = i;
    }
    ntbl_cnt = 0;
    dtbl_cnt = 0;
}

void update_closure() {
    for (int i = 0; i < ntbl_cnt; ++i) {
        int visited[TBL_SZ] = { 0 };
        MyQ<int> sq;
        sq.push(i);
        while (sq.size()) {
            int tmp = sq.front(); sq.pop();
            visited[tmp] = 1;
            ntbl[i].closure.insert(tmp);
            int n = ntbl[tmp].ne.size();
            for (int j = 0; j < n; ++j) {
                if (ntbl[tmp].ne[j].input == 'e' && visited[ntbl[tmp].ne[j].to->idx] == 0) {
                    sq.push(ntbl[tmp].ne[j].to->idx);
                }
            }
        }
    }
}

int is_new(IntSet& tgt) {
    for (int i = 0; i < dtbl_cnt; ++i) {
        if (tgt == dtbl[i].closure) return i;
    }
    return dtbl_cnt;
}

IntSet get_accessible_nodes(Node* ptr, char c) {
    IntSet ans;
    IntSet tmp;
    IntSet::iterator iter;
    for (iter = ptr->closure.begin(); iter != ptr->closure.end(); ++iter) {
        int n = ntbl[*iter].ne.size();
        for (int i = 0; i < n; ++i) {
            if (ntbl[*iter].ne[i].input == c) {
                tmp.insert(ntbl[*iter].ne[i].to->idx);
            }
        }
    }
    for (iter = tmp.begin(); iter != tmp.end(); ++iter) {
        ans.insert(*iter);
        ans.insert(ntbl[*iter].closure.begin(), ntbl[*iter].closure.end());
    }
    return ans;
}

Node* nfa2dfa(Graph& nfa) {
    MyQ<Node*> pq;
    Node* dptr = &dtbl[dtbl_cnt++];
    dptr->closure = nfa.start->closure;
    pq.push(dptr);
    Node* ans = dptr;
    MyVec<char> alphabet; alphabet.push_back('a'); alphabet.push_back('b');
    int n = alphabet.size();
    while (pq.size()) {
        Node* dptr = pq.front(); pq.pop();
        for (int i = 0; i < n; ++i) {
            IntSet tmp = get_accessible_nodes(dptr, alphabet[i]);
            if (tmp.size()) {
                int state_idx = is_new(tmp);
                if (state_idx == dtbl_cnt) {
                    Node* tmpptr = &dtbl[dtbl_cnt++];
                    tmpptr->closure = tmp;
                    pq.push(tmpptr);
                    Edge edge;
                    edge.input = alphabet[i];
                    edge.to = tmpptr;
                    dptr->ne.push_back(edge);
                }
                else {
                    Edge edge;
                    edge.input = alphabet[i];
                    edge.to = &dtbl[state_idx];
                    dptr->ne.push_back(edge);
                }
            }
        }
    }
    return ans;
}

void build_transition_mtx(Node* dfa, int **mtx) {
    MyQ<Node*> pq;
    pq.push(dfa);
    int visited[TBL_SZ] = { 0 };
    while (pq.size()) {
        Node* ptr = pq.front(); pq.pop();
        visited[ptr->idx] = 1;
        int n = ptr->ne.size();
        for (int i = 0; i < n; ++i) {
            mtx[ptr->idx][ptr->ne[i].to->idx] = 1;
            if (visited[ptr->ne[i].to->idx]==0) pq.push(ptr->ne[i].to);
        }
    }
}

void update_accepting(Graph &nfa) {
    for (int i = 0; i < dtbl_cnt; ++i) {
        if (dtbl[i].closure.find(nfa.end->idx) != dtbl[i].closure.end()) {
            dtbl[i].accept = true;
        }
    }
}

int count_strings(int dim, int** mod_mtx) {
    long long ans = 0;
    for (int i = 0; i < dim; ++i) {
        if (dtbl[i].accept) ans += mod_mtx[0][i];
    }
    return ans % MOD;
}

int countStrings(string &r, int l) {
    // initialize all global variables for reuse
    init_tbl();
    // convert the input regular expression from infix to postfix
    string re = in2post(r);
    // build a NFA from the postfix regular expression
    Graph nfa = build_nfa(re);
    // identify the closure for each state of the NFA
    update_closure();
    // build a DFA from the NFA
    Node* dfa = nfa2dfa(nfa);
    // identify the accepting states from the DFA
    update_accepting(nfa);
    // build a transition matrix
    int** mtx = alloc_sqr_mtx(dtbl_cnt);
    int** mod_mtx = alloc_sqr_mtx(dtbl_cnt);
    build_transition_mtx(dfa, mtx);
    // multiply the transition matrix l times
    mtx_pwr(dtbl_cnt, l, mtx, mod_mtx);
    // count all the transitions from the starting state with l characters
    int ans = count_strings(dtbl_cnt, mod_mtx);
    free_sqr_mtx(dtbl_cnt, mtx);
    free_sqr_mtx(dtbl_cnt, mod_mtx);
    return ans;
}

int main()
{
    ifstream ins("input.txt");
    ifstream res("result.txt");
    int TC; ins >> TC;
    for (int tc = 0; tc < TC; ++tc) {
        string reg_in; ins >> reg_in;
        int len; ins >> len;
        int ans = countStrings(reg_in, len);
        int truth; res >> truth;
        if (ans == truth) {
            cout << tc << "th test succeeded!" << endl;
        }
        else {
            cout << tc << "th test failed!" << endl;
            cout << "ans=" << ans << ", truth=" << truth << endl;
        }
    }
}
