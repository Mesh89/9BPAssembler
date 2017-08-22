/*
 * main.cpp
 *
 *  Created on: 10 Aug 2017
 *      Author: ramesh
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <stack>

using namespace std;

typedef unsigned long long ull;
typedef unsigned char uchar;

const int k = 30;
const double error_rate = 0.04;

ull nucl_bm[256] = { 0 }, rc_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };

int nucls_to_ulls(int nucls) {
	return nucls/32 + (nucls%32 > 0);
}

struct kmer_t {

	uchar nucls;
	ull* kmer;

	kmer_t() : nucls(0), kmer(NULL) {}

	kmer_t(string& s, bool rc = false) {
		nucls = s.size();
		int ints = nucls_to_ulls(nucls);
		kmer = new ull[ints];
		int offset = 62, intid = 0;
		kmer[0] = 0;
		for (int i = 0; i < nucls; i++) {
			if (offset == -2) {
				offset = 62;
				intid++;
				kmer[intid] = 0;
			}
			kmer[intid] |= ((rc ? rc_bm[s[nucls-1-i]] : nucl_bm[s[i]]) << offset);
			offset -= 2;
		}
	}

	char get_nucl(int pos) {
		return bm_nucl[((kmer[pos/32] >> (2*(32-pos%32)-2))) & 3];
	}

	string to_bin_str() {
		char* str = new char[2*nucls+1];
		ull mask = 1ULL << 63;
		for (int i = 0; i < nucls*2; i++) {
			str[i] = '0' + bool(kmer[i/64] & mask);
			mask >>= 1;
			if (mask == 0) mask = 1ULL << 63;
		}
		str[2*nucls] = '\0';
		return string(str);
	}

	string to_str() {
		char* str = new char[nucls+1];
		int shift = 62;
		for (int i = 0; i < nucls; i++) {
			str[i] = bm_nucl[(kmer[i/32] >> shift) & 3];
			shift -= 2;
		}
		str[nucls] = '\0';
		return string(str);
	}

	void shift_left_and_insert(char c) {
		// shift everything
		int last = nucls_to_ulls(nucls);
		for (int i = 0; i < last-1; i++) {
			kmer[i] <<= 2;
			kmer[i] |= kmer[i+1] >> 62;
		}
		kmer[last-1] <<= 2;

		// insert new nucl
		int zeroes = 2*(32-nucls%32);
		kmer[last-1] |= nucl_bm[c] << zeroes;
	}

	void shift_right_and_insert(char c) {
		int last = nucls_to_ulls(nucls)-1;

		int zeroes = 2*(32-nucls%32);
		kmer[last] >>= 2;
		kmer[last] = (kmer[last] >> zeroes) << zeroes;

		// shift everything
		for (int i = last-1; i >= 0; i--) {
			kmer[i+1] |= (kmer[i] & 3) << 62;
			kmer[i] >>= 2;
		}

		// insert new nucl
		kmer[0] |= nucl_bm[c] << 62;
	}

	kmer_t* prefix(int nucls) {
		kmer_t* prefix_kmer = new kmer_t;
		prefix_kmer->nucls = nucls;
		int ints = nucls_to_ulls(nucls);
		prefix_kmer->kmer = new ull[ints];
		memcpy(prefix_kmer->kmer, kmer, sizeof(ull)*ints);
		int zeroes = 2*(32-nucls%32);
		prefix_kmer->kmer[ints-1] = (prefix_kmer->kmer[ints-1] >> zeroes) << zeroes;
		return prefix_kmer;
	}

	~kmer_t() {
		delete[] kmer;
	}
};

struct kmer_p_hash {
	size_t operator()(const kmer_t* kmer) const {
		return hash<ull>() (kmer->kmer[0]);
	}
};

int kmer_p_cmp(kmer_t* const k1, kmer_t* const k2) {
	int nucls = k1->nucls < k2->nucls ? k1->nucls : k2->nucls;
	int ints = nucls_to_ulls(nucls);
	for (int i = 0; i < ints; i++) {
		if (k1->kmer[i] < k2->kmer[i]) return -1;
		else if (k1->kmer[i] > k2->kmer[i]) return 1;
	}
	return k1->nucls - k2->nucls;
}

struct kmer_p_eq {
	bool operator () (kmer_t* const k1, kmer_t* const k2) const {
		return kmer_p_cmp(k1, k2) == 0;
	}
};

bool kmer_p_less(kmer_t* const k1, kmer_t* const k2) {
	return kmer_p_cmp(k1, k2) < 0;
}

struct pair_t {
	kmer_t* r1,* r2;
	kmer_t* first1,* last1,* first2,* last2;

	pair_t(kmer_t* r1, kmer_t* r2) : r1(r1), r2(r2), first1(NULL), last1(NULL), first2(NULL), last2(NULL) {}
};

typedef unordered_map<kmer_t*, vector<pair<pair_t*, uchar> >, kmer_p_hash, kmer_p_eq> kmer_indexed_pairs;
typedef unordered_map<kmer_t*, vector<pair<pair_t*, pair<uchar, uchar> > >, kmer_p_hash, kmer_p_eq> kmer_indexed_pairs_2pos;

kmer_indexed_pairs unique_kmers;


bool contains_N(string& r) {
	return r.find_first_of('N') != string::npos;
}

void index_read1(string& r, pair_t* p) {
	string temp = r.substr(0, k);
	kmer_t* ukmer = new kmer_t(temp);
	uchar pos = 0;
	for (int i = k; i < r.length(); i++) {
		auto it = unique_kmers.find(ukmer);
		if (it != unique_kmers.end()) {
			if (p->first1 == NULL) p->first1 = it->first;
			p->last1 = it->first;
			it->second.push_back(make_pair(p, pos));
		}
		if (i == r.length()) {
			ukmer = NULL;
			break;
		}
		ukmer->shift_left_and_insert(r[i]);
		pos++;
	}
	delete ukmer;
}

bool kmer_pos_less(const pair<kmer_t*, uchar>& p1, const pair<kmer_t*, uchar>& p2) {
	return kmer_p_less(p1.first, p2.first);
}

// TODO: avoid read -> string -> kmers by using kmer->prefix
// Given a read, returns its unique kmers along with their position
void index_read2(pair_t* p, uchar pos1, kmer_indexed_pairs_2pos& r2index) {
	kmer_t* read = p->r2;
	string read_str = read->to_str();
	string temp = read_str.substr(0, k);
	kmer_t* ukmer = new kmer_t(temp);
	uchar pos = 0;
	for (int i = k; i < read->nucls; i++) {
		auto it = unique_kmers.find(ukmer);
		if (it != unique_kmers.end()) {
			if (p->first2 == NULL) p->first2 = it->first;
			p->last2 = it->first;
			r2index[it->first].push_back(make_pair(p, make_pair(pos1, pos)));
		}
		ukmer->shift_left_and_insert(read_str[i]);
		pos++;
	}
}

double get_score(kmer_t* r1, uchar pos1, kmer_t* r2, uchar pos2) {
	int matches = 0, overlap = 0;
	while (pos1 < r1->nucls && pos2 < r2->nucls) {
		matches += r1->get_nucl(pos1++) == r2->get_nucl(pos2++);
		overlap++;
	}
	return double(matches)/overlap;
}

double get_score(pair_t* p1, pair<uchar, uchar> p1_pos, pair_t* p2, pair<uchar, uchar> p2_pos) {
	return get_score(p1->r1, p1_pos.first, p2->r1, p2_pos.first) * get_score(p1->r2, p1_pos.second, p2->r2, p2_pos.second);
}

void get_component(pair_t* p, unordered_map<pair_t*, vector<pair_t*> >& edges, unordered_set<pair_t*>& pairs_used, vector<pair_t*>& comp) {
	stack<pair_t*> s;
	s.push(p);
	pairs_used.insert(p);
	while (!s.empty()) {
		p = s.top();
		comp.push_back(p);
		s.pop();
		for (pair_t* dest : edges[p]) {
			if (!pairs_used.count(dest)) {
				s.push(dest);
				pairs_used.insert(dest);
			}
		}
	}
}

int main(int argc, char* argv[]) {

	if (argc != 4) {
		cout << "Usage: ./exec unique-kmers file1.fa file2.fa" << endl;
		return 0;
	}

	ifstream ukmersfile(argv[1]), fa1(argv[2]), fa2(argv[3]);

	nucl_bm['A'] = 0, nucl_bm['C'] = 1, nucl_bm['G'] = 2, nucl_bm['T'] = 3;
	rc_bm['T'] = 0, rc_bm['G'] = 1, rc_bm['C'] = 2, rc_bm['A'] = 3;


	string ukmer;
	while (getline(ukmersfile, ukmer)) {
		if (ukmer.length() != k) {
			cout << "k-mer length not " << k << endl;
			return -1;
		}
		kmer_t* kmer = new kmer_t(ukmer);
		unique_kmers[kmer] = vector<pair<pair_t*, uchar> >();

		kmer = new kmer_t(ukmer, true);
		unique_kmers[kmer] = vector<pair<pair_t*, uchar> >();
	}


	unordered_map<pair_t*, string> pair_ids; // TODO: remove

	string id1, r1, id2, r2;
	vector<pair_t*> read_pairs;
	while (getline(fa1, id1) && getline(fa2, id2)) {
		getline(fa1, r1); getline(fa2, r2);

		if (contains_N(r1) || contains_N(r2)) continue;

		pair_t* p = new pair_t(new kmer_t(r1), new kmer_t(r2));
		read_pairs.push_back(p);
		index_read1(r1, p);
		pair_ids[p] = id1.substr(1, id1.length()-3);
	}

	unordered_map<pair_t*, vector<pair_t*> > edges;

	int seen = 0;
	for (pair<kmer_t*, vector<pair<pair_t*, uchar> > > pairs_with_ukmer1 : unique_kmers) {
		if (pairs_with_ukmer1.second.empty()) continue;

		kmer_indexed_pairs_2pos index_by_r2;
		// index reads sharing a unique k-mer in reads 1 by unique k-mers in reads 2
		for (pair<pair_t*, uchar> pair_with_ukmer1 : pairs_with_ukmer1.second) {
			pair_t* p = pair_with_ukmer1.first;
			index_read2(p, pair_with_ukmer1.second, index_by_r2);
		}

		// process pairs sharing a unique k-mer in both reads
		for (pair<kmer_t*, vector<pair<pair_t*, pair<uchar, uchar> > > > _pairs_with_ukmer2 : index_by_r2) {
			vector<pair<pair_t*, pair<uchar, uchar> > >& pairs_with_ukmer2 = _pairs_with_ukmer2.second;
			for (int i = 0; i < pairs_with_ukmer2.size(); i++) {
				pair_t* pair1 = pairs_with_ukmer2[i].first;
				if ((pair1->first1 != pairs_with_ukmer1.first && pair1->last1 != pairs_with_ukmer1.first) ||
					(pair1->first2 != _pairs_with_ukmer2.first && pair1->last2 != _pairs_with_ukmer2.first)) continue;
				pair<uchar, uchar> pair1_pos = pairs_with_ukmer2[i].second;
				for (int j = i+1; j < pairs_with_ukmer2.size(); j++) {
					pair_t* pair2 = pairs_with_ukmer2[j].first;
					pair<uchar, uchar> pair2_pos = pairs_with_ukmer2[j].second;
					double score = get_score(pair1, pair1_pos, pair2, pair2_pos);
					if (score >= 1-error_rate) {
						edges[pair1].push_back(pair2);
					}
				}
			}
		}

		if (seen % 10000 == 0) cerr << double(seen)/unique_kmers.size()*100 << "%" << endl;
		seen++;
	}

	unordered_set<pair_t*> pairs_used;

	for (pair_t* p : read_pairs) {
		if (pairs_used.count(p)) continue;
		vector<pair_t*> comp;
		get_component(p, edges, pairs_used, comp);
		cout << comp.size() << endl;
	}
}
