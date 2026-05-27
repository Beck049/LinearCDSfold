#ifndef _AC_MATCHER_H_
#define _AC_MATCHER_H_

#include <cstring>
#include <list>
#include <queue>
#include <unordered_map>

namespace oligodesign {

/**
 * @brief Aho-Corasick 字串匹配演算法的實作，支援多模式字串搜尋。
 * 
 * @tparam T 字元型態（例如 char）
 * @tparam S 字元集大小（例如 256）
 */
template<class T, int S> 
class ACMatcher {
public:
    ACMatcher();
    ~ACMatcher();
    void AddString(const T *, int);
    bool MakeTree();
    int Search(const T *, int);
    int Search(const T *, int, int);
    void SearchPositions(const T *text,int idx,int start,int end, std::unordered_map<int,int> &mp);
    int NumEndingInLast3(const T *, int);
    int NumEndingInLastK(const T *, int, int, int);

private:
    struct ACState {
        int id;
        int depth;
        int count;
        T *output;
        struct ACState *fail;
        struct ACState *transitions[S];
    };
    typedef struct ACState ACState_t;
      
    int newstate;
    ACState_t *zerostate;
    std::list<T *> cached_strings;
        
    void ACFree(ACState_t *);
};

/**
 * @brief 建構子：初始化自動機
 */
template<class T, int S>
ACMatcher<T,S>::ACMatcher() {
    zerostate = 0;
    newstate = 0;
}

/**
 * @brief 解構子：釋放樹狀結構與所有快取字串的記憶體
 */
template<class T, int S>
ACMatcher<T,S>::~ACMatcher() {
    if (zerostate != 0) {
        for (int i = 0; i < S ;i++) {
            if (zerostate->transitions[i] != 0 &&
                zerostate->transitions[i]->id > 0) {
                ACFree(zerostate->transitions[i]);
            }
        }
    }
    // clean up cached strings
    typename std::list<T *>::iterator itr;
    for( itr = cached_strings.begin(); itr != cached_strings.end(); itr++ ) {
        if (*itr != 0) {
            delete [] *itr;
        }
    }
}

/**
 * @brief 輔助函式：遞迴釋放特定節點及其子樹的記憶體
 */
template<class T, int S>
void ACMatcher<T,S>::ACFree(ACState_t *state) {
    for (int i = 0; i < S ;i++) {
        if (state->transitions[i] != 0) {
            ACFree(state->transitions[i]);
        }
    }
    delete(state);
}
/**
 * @brief 新增一個關鍵字（模式串）到字典樹（Trie Tree）中
 * @param text 關鍵字陣列指標
 * @param n 關鍵字長度
 */
template<class T, int S>
void ACMatcher<T,S>::AddString(const T *text, int n) {
    if (text == 0) return;
      
    // first cache the string locally
    cached_strings.push_back(new T[n]);
    T *cached_item = cached_strings.back();
    for(int i=0; i<n; i++) {
        cached_item[i] = text[i];
    }
      
    ACState_t *state = 0, *s = 0;
    int j = 0;

    if (zerostate == 0) {
        zerostate = new ACState_t;
        newstate = 1;
        zerostate->id = 0;
        zerostate->depth = 0;
        zerostate->output = 0;
        memset(zerostate->transitions, 0, S * sizeof(ACState_t *));
    }

    state = zerostate;

    // As long as we have transitions follow them
    while (j < n && (s = state->transitions[(int)*(text+j)]) != 0) {
        state = s;
        ++j;
    }

    if (j == n) {
        s->count = s->count + 1;
        return;
    }

    while (j < n) {
        // Create new state
        s = new ACState_t;
        s->id = newstate++;
        s->depth = state->depth + 1;
        memset(s->transitions, 0, S * sizeof(ACState_t *));
        // Create transition
        state->transitions[(int) * (text + j)] = s;
        state = s;
        s->output = 0;
        s->count = 1;
        ++j;
    }

    s->output = cached_strings.back();
    return;
}

/**
 * @brief 建立失敗指標（Fail Pointers）。必須在 AddString 全部完成後呼叫。
 * @note 透過 BFS (廣度優先搜尋) 來串接每個節點的 fail 指標
 * @return bool 建立成功傳回 true，若無字串則傳回 false
 */
template<class T, int S>
bool ACMatcher<T,S>::MakeTree() {
    if (cached_strings.empty()) return false;
  
    std::queue<ACState_t *> state_queue;
    ACState_t *state, *s, *r;
    int i;

    // Set all FAIL transition of 0 state to point to itself
    for (i = 0; i < S; i++) {
        if (zerostate->transitions[i] == 0)
            zerostate->transitions[i] = zerostate;
            // Construct fail()
        else {
            state_queue.push(zerostate->transitions[i]);
            zerostate->transitions[i]->fail = zerostate;
        }
    }

    // Set fail() for depth > 0
    while (!state_queue.empty()) {
        r = state_queue.front();
        state_queue.pop();
        for (i = 0; i < S; i++) {
            if ((s = r->transitions[i]) == 0) {
                continue;
            }
            state_queue.push(s);
            state = r->fail;
            while (state->transitions[i] == 0) {
                state = state->fail;
            }
            s->fail = state->transitions[i];
            // Join outputs missing
        }
    }
    return true;
}

/**
 * @brief 在整個文本中搜尋關鍵字
 * @param text 待搜尋的文本
 * @param n 文本長度
 * @return int 總匹配成功的次數
 */
template<class T, int S>
int ACMatcher<T,S>::Search(const T *text, int n) {
    if (text == 0) return 0;
    return Search(text, 0, n - 1);
}

/**
 * @brief 在文本的指定範圍內搜尋關鍵字 [start, end]
 * @param text 待搜尋的文本
 * @param start 起始索引
 * @param end 結束索引
 * @return int 總匹配成功的次數
 */
template<class T, int S>
int ACMatcher<T,S>::Search(const T *text, int start, int end) {
    if (cached_strings.empty()) return 0;
  
    ACState_t *state = zerostate;
    int num_matches = 0;

    for (int j = start; j <= end; j++) {
        while (state->transitions[(int)*(text + j)] == 0) {
            state = state->fail;
        }
        state = state->transitions[(int)*(text + j)];
        if (state->output != 0) {
            //return state->id;
            num_matches += state->count;
            //return j;
        }
    }
    return num_matches;
    //return -1;
}

/**
 * @brief 搜尋並記錄匹配到關鍵字的「結束位置」與「長度」
 * @param text 待搜尋的文本
 * @param idx 偏移量（用於計算全域位置）
 * @param start 搜尋起點
 * @param end 搜尋終點
 * @param mp 用於儲存結果的 map，Key 為文字結束位置 (j + idx)，Value 為匹配字串的長度 (depth)
 */
template<class T, int S>
void ACMatcher<T,S>::SearchPositions(const T *text,int idx,int start,int end, std::unordered_map<int,int> &mp) {
    if (cached_strings.empty() || text == 0)
        return ;
    ACState_t *state = zerostate;
    for (int j = start; j <= end; j++) {
        while (state->transitions[(int)*(text + j)] == 0) {
            state = state->fail;
        }
        state = state->transitions[(int)*(text + j)];
        if (state->output != 0) {
            for (int k = 0; k < state->count; k++) {
                if(mp.find(idx+j) == mp.end()){
                    mp[j+idx] = state->depth;
                }
            }
        }
    }
    return;
}

/**
 * @brief 統計在文本「最後 3 個字元內」結尾的關鍵字數量
 */
template<class T, int S>
int ACMatcher<T,S>::NumEndingInLast3(const T *text, int n) {
    if (text == 0) return 0;
    return NumEndingInLastK(text, 0, n - 1, 3);
}

/**
 * @brief 統計在文本「最後 K 個字元內」結尾的關鍵字數量
 * @param text 待搜尋的文本
 * @param start 搜尋起點
 * @param end 搜尋終點
 * @param k 範圍限制（只統計結束索引大於 end - k 的匹配）
 * @return int 符合條件的匹配數量
 */
template<class T, int S>
int ACMatcher<T,S>::NumEndingInLastK(const T *text, int start, int end, int k) {
    if (cached_strings.empty()) return 0;
    ACState_t *state = zerostate;
    int num_matches = 0;

    for (int j = start; j <= end; j++) {
        while (state->transitions[(int)*(text + j)] == 0) {
            state = state->fail;
        }
        state = state->transitions[(int)*(text + j)];
        if (state->output != 0) {
            //return state->id;
            if (j > end - k) {
                num_matches += state->count;
            }
            //return j;
        }
    }
    return num_matches;
    //return -1;
}

} // end namespace oligodesign

#endif /* _AC_MATCHER_H_ */
