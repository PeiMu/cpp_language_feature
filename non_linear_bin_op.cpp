#include <algorithm>
#include <assert.h>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

bool clearly = false;

std::pair<int64_t, int64_t>
modInterval(const int64_t lhs_low, const int64_t lhs_high, const int64_t rhs_low, const int64_t rhs_high);

std::pair<int64_t, int64_t>
modConstDivisorInterval(const int64_t lhs_low, const int64_t lhs_high, const int64_t m)
{
	// empty interval
	if (lhs_low > lhs_high || m == 0)
		return std::make_pair(0, 0);
	// compute modulo with positive interval and negate
	else if (lhs_high < 0)
		return std::make_pair(-modConstDivisorInterval(-lhs_high, -lhs_low, m).second,
                              -modConstDivisorInterval(-lhs_high, -lhs_low, m).first);
	// split into negative and non-negative interval, compute and join
	else if (lhs_low < 0)
	{
		auto negative_part = modConstDivisorInterval(lhs_low, -1, m);
		auto positive_part = modConstDivisorInterval(0, lhs_high, m);
		return std::make_pair(std::min(negative_part.first, positive_part.first),
				      std::max(negative_part.second, positive_part.second));
	}
	// there is no k > 0 such that a < k*m <= b
	else if ((lhs_high - lhs_low) < std::abs(m) && (lhs_low % m <= lhs_high % m))
		return std::make_pair(lhs_low % m, lhs_high % m);
	else
		return std::make_pair(0, std::abs(m) - 1);
}

std::pair<int64_t, int64_t>
modConstDividendInterval(const int64_t lhs, const int64_t rhs_low, const int64_t rhs_high)
{
    // empty interval
    if (rhs_low > rhs_high || lhs == 0)
        return std::make_pair(0, 0);
    // compute modulo with positive interval and negate
    else if (lhs < 0)
        return std::make_pair(-modConstDividendInterval(-lhs, rhs_low, rhs_high).second,
                              -modConstDividendInterval(-lhs, rhs_low, rhs_high).first);
    // use only non-negative m and n
    else if (rhs_high <= 0)
        return modConstDividendInterval(lhs, -rhs_high, -rhs_low);
    // split into negative and non-negative interval, compute and join
    else if (rhs_low <= 0)
        return modConstDividendInterval(lhs, 1, std::max(-rhs_low, rhs_high));
    // modulo has no effect
    else if (rhs_low > lhs)
        return std::make_pair(lhs, lhs);
    // there is some overlapping of [a,b] and [n,m]
    else if (rhs_high > lhs)
        return std::make_pair(0, lhs);
    // max value appears at [a/2] + 1
    else if (lhs / 2 + 1 <= rhs_low)
        return std::make_pair(lhs % rhs_high, lhs % rhs_low);
    else if (lhs / 2 + 1 <= rhs_high) {
        int64_t min = lhs;
        if (rhs_high != lhs) {
            for (int64_t min_idx = rhs_low; min_idx < lhs / 2 + 1; min_idx++) {
                min = lhs % min_idx < min ? lhs % min_idx : min;
            }
        } else {
            min = 0;
        }
        return std::make_pair(min, lhs % (lhs / 2 + 1));
    }
    // either compute all possibilities and join, or be imprecise
    else
    {
        int64_t min = lhs;
        int64_t max = 0;
        for (int64_t min_idx = rhs_low; min_idx <= rhs_high; min_idx++) {
            min = lhs % min_idx < min ? lhs % min_idx : min;
        }
        for (int64_t max_idx = rhs_low; max_idx <= rhs_high; max_idx++) {
            max = lhs % max_idx > max ? lhs % max_idx : max;
        }
        return std::make_pair(min, max);
    }
}

std::pair<int64_t, int64_t>
modInterval(const int64_t lhs_low, const int64_t lhs_high, const int64_t rhs_low, const int64_t rhs_high)
{
    // empty interval
    if (lhs_low > lhs_high || rhs_low > rhs_high)
        return std::make_pair(0, 0);
    // compute modulo with positive interval and negate
    else if (lhs_high < 0)
        return std::make_pair(-modInterval(-lhs_high, -lhs_low, rhs_low, rhs_high).second,
                              -modInterval(-lhs_high, -lhs_low, rhs_low, rhs_high).first);
    // split into negative and non-negative interval, compute, and join
    else if (lhs_low < 0)
    {
        auto negative_part = modInterval(lhs_low, -1, rhs_low, rhs_high);
        auto positive_part = modInterval(0, lhs_high, rhs_low, rhs_high);
        return std::make_pair(std::min(negative_part.first, positive_part.first),
                              std::max(negative_part.second, positive_part.second));
    }
    // use the simpler function from before
    else if (lhs_low == lhs_high)
        return modConstDividendInterval(lhs_low, rhs_low, rhs_high);
    // use the simpler function from before
    else if (rhs_low == rhs_high)
        return modConstDivisorInterval(lhs_low, lhs_high, rhs_low);
    // use only non-negative m and n
    else if (rhs_high <= 0)
        return modInterval(lhs_low, lhs_high, -rhs_high, -rhs_low);
    // make modulus non-negative
    else if (rhs_low <= 0)
        return modInterval(lhs_low, lhs_high, 1, std::max(-rhs_low, rhs_high));
    // check b-a < |modulus|
    else if (lhs_high - lhs_low >= rhs_high)
        return std::make_pair(0, rhs_high - 1);
    // split interval, compute, and join
    else if (lhs_high - lhs_low >= rhs_low)
    {
        auto part = modInterval(lhs_low, lhs_high, lhs_high-lhs_low+1, rhs_high);
        return std::make_pair(std::min((int64_t)0, part.first), std::max(lhs_high-lhs_low-1, part.second));
    }
    // modulo has no effect
    else if (rhs_low > lhs_high)
        return std::make_pair(lhs_low, lhs_high);
    // there is some overlapping of [a,b] and [n,m]
    else if (rhs_high > lhs_high)
        return std::make_pair(0, lhs_high);
    // either compute all possibilities and join, or be imprecise
    else {
        auto dist_lhs = lhs_high - lhs_low + 1;
        auto dist_rhs = rhs_high - rhs_low + 1;
        std::vector<std::pair<int64_t, int64_t>> res;
        if (dist_lhs < dist_rhs) {
            for (int64_t lhs = lhs_low; lhs <= lhs_high; lhs++) {
                res.emplace_back(modConstDividendInterval(lhs, rhs_low, rhs_high));
            }
        } else {
            for (int64_t rhs = rhs_low; rhs <= rhs_high; rhs++) {
                res.emplace_back(modConstDivisorInterval(lhs_low, lhs_high, rhs));
            }
        }
        auto min_res = std::min_element(res.begin(), res.end(), [](auto a, auto b) {
            return a.first < b.first;
        });
        auto max_res = std::max_element(res.begin(), res.end(), [](auto a, auto b) {
            return a.second < b.second;
        });
        return std::make_pair(min_res->first, max_res->second);
    }

}

void
checkModRes(const int64_t lhs_low, const int64_t lhs_high,
            const int64_t rhs_low, const int64_t rhs_high,
            std::ofstream& ofs)
{
    if (rhs_low == 0)
        return;
	int64_t rem_low	 = lhs_low % rhs_low;
	int64_t rem_high = lhs_low % rhs_low;
    #pragma omp parallel for
	for (int64_t i = lhs_low; i <= lhs_high; i++)
	{
        #pragma omp parallel for
		for (int64_t j = rhs_low; j <= rhs_high; j++)
		{
			if (j == 0)
                continue;
			int64_t tmp_rem = i % j;
            if (clearly)
                std::cout << "tmp_rem: " << i << " % " << j << " = " << tmp_rem << std::endl;
			rem_low	 = tmp_rem < rem_low ? tmp_rem : rem_low;
			rem_high = tmp_rem > rem_high ? tmp_rem : rem_high;
		}
	}
    std::pair<int64_t, int64_t> compute_res = modInterval(lhs_low, lhs_high, rhs_low, rhs_high);
    if (compute_res.first == rem_low && compute_res.second == rem_high) {
//        printf("compute success with[%ld, %ld] mod [%ld, %ld]. real range is [%ld, %ld]\n",
//               lhs_low, lhs_high, rhs_low, rhs_high, rem_low, rem_high);
        ofs << "compute success with [" << lhs_low << ", " << lhs_high << "] % [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << rem_low << ", " << rem_high << "]" << std::endl;
    }
    else
    {
//        printf("compute error with [%ld, %ld] mod [%ld, %ld]. real range is [%ld, %ld], "
//               "compute range is [%ld, %ld]\n", lhs_low, lhs_high, rhs_low, rhs_high, rem_low,
//               rem_high, compute_res.first, compute_res.second);
        ofs << "compute error with [" << lhs_low << ", " << lhs_high << "] % [" << rhs_low << ", "
        << rhs_high << "]. real range is [" << rem_low << ", " << rem_high << "], compute range is ["
        << compute_res.first << ", " << compute_res.second << "]" << std::endl;
    }
}

std::vector<std::string> getDisjointUnionSet(const boost::dynamic_bitset<>& lower_bin_set,
                                             const boost::dynamic_bitset<>& upper_bin_set) {
    std::string lower_bin, upper_bin;
    to_string(lower_bin_set, lower_bin);
    to_string(upper_bin_set, upper_bin);
    const uint8_t bit_size = lower_bin.length();
    std::vector<std::string> disjoint_union_set;
    std::string bin_str_ele;

    /*
     * the turning point is where lower_bin different from upper_bin
     * */
    std::string xor_str;
    to_string((lower_bin_set ^ upper_bin_set), xor_str);
    if (xor_str.find_first_of('1') == std::string::npos) {
        /*
         * lower bound == higher bound, return with itself
         * */
        disjoint_union_set.emplace_back(lower_bin);
        return disjoint_union_set;
    }
    const auto turning_point = xor_str.find_first_of('1');
    /*
     * for the right-most continuously bit of lower bound,
     *  if they're '0', then change them to 'x'
     *  if it's '1', then emplace '1'
     */
    const auto lower_start_pos = lower_bin.substr(1).find_last_of('1') + 1;
    if (lower_start_pos == bit_size-1) {
        disjoint_union_set.emplace_back(lower_bin);
    } else {
        /* substitute all of the right most '0' in lower_bin to 'x'
         *  e.g.
         *  ([010100], [011101]) -> {[0101xx], ...}
         * but if the turning point is righter than 'lower_start_pos'
         * only substitute bits right after the turning point
         *  e.g.
         *  ([01000], [01110]) -> {[010xx], ...}
         * */
        const auto substitute_bit = lower_start_pos > turning_point ? lower_start_pos : turning_point;
        bin_str_ele = lower_bin.substr(0, substitute_bit + 1);
        bin_str_ele.append(bit_size-substitute_bit - 1, 'x');
        disjoint_union_set.emplace_back(bin_str_ele);
        bin_str_ele.clear();
    }
    /*
     * when splitting the lower bound, iterate from right to left
     * */
    for (int idx = lower_start_pos-1; idx > turning_point; idx--) {
        if (idx < 0)
            break;
        /*
         * substitute the '0' to '1' and append with 'x' to get the large number
         * */
        if (lower_bin[idx] == '0') {
            bin_str_ele = lower_bin.substr(0, idx) + '1';
            bin_str_ele.append(bit_size-idx-1, 'x');
            disjoint_union_set.emplace_back(bin_str_ele);
            bin_str_ele.clear();
        }
    }

    /*
     * for the right-most continuously bit of upper bound,
     *  if they're '1', then change them to 'x'
     *  if it's '0', then emplace '0'
     * */
    const auto upper_end_pos = upper_bin.substr(1).find_last_of('0') + 1;
    if (upper_end_pos == bit_size-1) {
        disjoint_union_set.emplace_back(upper_bin);
    } else {
        /* substitute all of the right most '1' in upper_bin to 'x'.
         *  e.g.
         *  ([010001], [011011]) -> {[0110xx], ...}
         * but if the turning point is righter than 'upper_end_pos'
         * only substitute bits right after the turning point
         *  e.g.
         *  ([010001], [011111]) -> {[011xxx], ...}
         * */
        const auto substitute_bit = upper_end_pos > turning_point ? upper_end_pos : turning_point;
        bin_str_ele = upper_bin.substr(0, substitute_bit + 1);
        bin_str_ele.append(bit_size-substitute_bit - 1, 'x');
        disjoint_union_set.emplace_back(bin_str_ele);
        bin_str_ele.clear();
    }
    /*
     * when splitting the upper bound, iterate from left to right
     * */
    for (int idx = turning_point+1; idx < upper_end_pos; idx++) {
        /*
         * substitute the '1' to '0' and append with 'x' to get the smaller number
         * */
        if (upper_bin[idx] == '1') {
            bin_str_ele = upper_bin.substr(0, idx) + '0';
            bin_str_ele.append(bit_size-idx-1, 'x');
            disjoint_union_set.emplace_back(bin_str_ele);
            bin_str_ele.clear();
        }
    }

    return disjoint_union_set;
}

enum BitSetType {
    HAS_BIT = 1,
    ALL_SAME = 2,
    ALL_DIFF = 3,
};

using bitSets = std::vector<std::vector<size_t>>;
using signedBit = std::pair<bool, std::vector<size_t>>;
struct hashFunc {
    size_t operator()(const signedBit& x) const {
        size_t second_xor;
        for (const auto& s : x.second) {
            if (s == SIZE_MAX)
                return x.first;
            second_xor ^= s;
        }
        return x.first ^ second_xor;
    }
};
using signedBitSets = std::unordered_set<signedBit, hashFunc>;

template <BitSetType BST, typename Op>
signedBitSets getBitPos(bitSets lhs_vecs, bitSets rhs_vecs,
                        std::vector<std::string> lhs_dus,
                        std::vector<std::string> rhs_dus,
                        Op bitwise_func) {
    signedBitSets res;
    const std::vector<size_t> non_exit(1, SIZE_MAX);
    for (size_t lhs_idx = 0; lhs_idx < lhs_vecs.size(); lhs_idx++) {
        for (size_t rhs_idx = 0; rhs_idx < rhs_vecs.size(); rhs_idx++) {
            std::vector<size_t> tmp_same_pos;
            if constexpr (BST == BitSetType::HAS_BIT) {
                std::set_union(lhs_vecs[lhs_idx].begin(), lhs_vecs[lhs_idx].end(),
                               rhs_vecs[rhs_idx].begin(), rhs_vecs[rhs_idx].end(),
                               std::back_inserter(tmp_same_pos));
            } else if constexpr (BST == BitSetType::ALL_SAME) {
                std::set_intersection(lhs_vecs[lhs_idx].begin(), lhs_vecs[lhs_idx].end(),
                                      rhs_vecs[rhs_idx].begin(), rhs_vecs[rhs_idx].end(),
                                      std::back_inserter(tmp_same_pos));
            } else if constexpr (BST == BitSetType::ALL_DIFF) {
                std::set_difference(lhs_vecs[lhs_idx].begin(), lhs_vecs[lhs_idx].end(),
                                    rhs_vecs[rhs_idx].begin(), rhs_vecs[rhs_idx].end(),
                                    std::back_inserter(tmp_same_pos));
            } else
                assert(false && "unknown enum of BitSetType");
            if (!tmp_same_pos.empty()) {
                res.emplace(std::make_pair(bitwise_func(lhs_dus[lhs_idx][0]-48,
                                                                      rhs_dus[rhs_idx][0]-48),
                                                         tmp_same_pos));
            } else {
                res.emplace(std::make_pair(bitwise_func(lhs_dus[lhs_idx][0]-48,
                                                                      rhs_dus[rhs_idx][0]-48),
                                                         non_exit));
            }
        }
    }
    return res;
}

int64_t convertBin2Dec(const std::string& min_res_str, const uint8_t bit_size) {
    int64_t res;
    if (min_res_str[0] == '1') {
        auto tmp_res = min_res_str;
        /*
         * minus 1
         * */
        if (tmp_res.back() == '1') {
            tmp_res.back() = '0';
        } else {
            const size_t last_one = tmp_res.find_last_of('1');
            std::string append_one(bit_size-last_one-1, '1');
            tmp_res[last_one] = '0';
            tmp_res = tmp_res.substr(0, last_one+1) + append_one;
        }
        /*
         * negate
         * */
        boost::dynamic_bitset<> min_res_db(tmp_res);
        res = -min_res_db.flip().to_ulong();
    } else {
        boost::dynamic_bitset<> min_res_db(min_res_str);
        res = min_res_db.to_ulong();
    }
    return res;
}

template <typename Op>
std::pair<int64_t, int64_t>
bitwiseInterval(const int64_t lhs_low, const int64_t lhs_high,
                const int64_t rhs_low, const int64_t rhs_high,
                Op bitwise_func)
{
    /*
     * 1. get the max bit size
     *  for negative number, it has the same or smaller remaining bit size with its positive,
     *  like -8="1...1000", 8="0...1000"; -14="1...0010", 14="0...1110"
     *  so when meet negative number, we get its max bit size with its absolute value.
     *  PS: the first bit is a sign bit.
     * */
    std::vector<int64_t> ranges{std::abs(lhs_low), std::abs(lhs_high),
                                std::abs(rhs_low), std::abs(rhs_high)};
    auto max_abs_value = std::max_element(ranges.begin(), ranges.end());
    uint8_t bit_size = 64 - std::bitset<64>(*max_abs_value).to_string().find('1') + 1;

    /*
     * 2. convert to binary number
     * */
    boost::dynamic_bitset<> ll_bin(bit_size, lhs_low);
    boost::dynamic_bitset<> lh_bin(bit_size, lhs_high);
    boost::dynamic_bitset<> rl_bin(bit_size, rhs_low);
    boost::dynamic_bitset<> rh_bin(bit_size, rhs_high);

    /*
     * 3. construct disjoint union set of each binary number
     * */
    auto lhs_dus = getDisjointUnionSet(ll_bin, lh_bin);
    auto rhs_dus = getDisjointUnionSet(rl_bin, rh_bin);

    /*
     * 4. for different bit_wise operation, we have different formulas and purpose
     * */
    /*collect all bits index with '0' and '1'*/
    bitSets lhs_zero_vecs, lhs_one_vecs;
    for (const auto& lhs: lhs_dus) {
        std::vector<size_t> match_zero, match_one;
        /*no need to check the sign bit*/
        std::vector<size_t> index_vec(bit_size-1);
        std::iota(index_vec.begin(), index_vec.end(), 1);
        std::copy_if(index_vec.begin(), index_vec.end(), std::back_inserter(match_zero), [lhs](size_t v){
            return lhs[v] == '0';
        });
        lhs_zero_vecs.emplace_back(match_zero);
        std::copy_if(index_vec.begin(), index_vec.end(), std::back_inserter(match_one), [lhs](size_t v){
            return lhs[v] == '1';
        });
        lhs_one_vecs.emplace_back(match_one);
    }
    bitSets rhs_zero_vecs, rhs_one_vecs;
    for (const auto& rhs: rhs_dus) {
        std::vector<size_t> match_zero, match_one;
        /*no need to check the sign bit*/
        std::vector<size_t> index_vec(bit_size-1);
        std::iota(index_vec.begin(), index_vec.end(), 1);
        std::copy_if(index_vec.begin(), index_vec.end(), std::back_inserter(match_zero), [rhs](size_t v){
            return rhs[v] == '0';
        });
        rhs_zero_vecs.emplace_back(match_zero);
        std::copy_if(index_vec.begin(), index_vec.end(), std::back_inserter(match_one), [rhs](size_t v){
            return rhs[v] == '1';
        });
        rhs_one_vecs.emplace_back(match_one);
    }

    signedBitSets same_one_pos = getBitPos<BitSetType::ALL_SAME>
            (lhs_one_vecs, rhs_one_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets same_zero_pos = getBitPos<BitSetType::ALL_SAME>
            (lhs_zero_vecs, rhs_zero_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets exist_one_pos = getBitPos<BitSetType::HAS_BIT>
            (lhs_one_vecs, rhs_one_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets exist_zero_pos = getBitPos<BitSetType::HAS_BIT>
            (lhs_zero_vecs, rhs_zero_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets diff_first_pos = getBitPos<BitSetType::ALL_SAME>
            (lhs_one_vecs, rhs_zero_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets diff_second_pos = getBitPos<BitSetType::ALL_SAME>
            (lhs_zero_vecs, rhs_one_vecs, lhs_dus, rhs_dus, bitwise_func);
    signedBitSets union_min, union_max;
    // todo: there's a bug of xor: [-5, -4, 0, 1]
    std::set_union(diff_first_pos.begin(), diff_first_pos.end(),
                   diff_second_pos.begin(), diff_second_pos.end(),
                   std::inserter(union_min, union_min.begin()));
    std::set_union(same_zero_pos.begin(), same_zero_pos.end(),
                   same_one_pos.begin(), same_one_pos.end(),
                   std::inserter(union_max, union_max.begin()));

    signedBitSets::iterator min_res_it, max_res_it;
    auto compSingBit = [](const signedBit& a, const signedBit& b) {
        auto it_a = a.second.begin();
        auto it_b = b.second.begin();
        while (it_a != a.second.end() && it_b != b.second.end()) {
            if ((*it_a) < (*it_b))
                return true;
            else if ((*it_a) > (*it_b)){
                return false;
            } else {
                it_a++;
                it_b++;
            }
        }
        return it_a!=a.second.end();
    };
    if constexpr (std::is_same< decltype(bitwise_func), decltype(std::bit_and<bool>())>::value) {
        /*
         * and:
         *  0 & 0 = 0, 0 & 1 = 0, 0 & x = 0
         *  1 & 1 = 1, 1 & x = x, x & x = x
         *   min: find the same '1', the larger the better, others are '0'
         *   max: find '0' exist, the larger the better, others are '1'
         * */
        min_res_it = std::max_element(same_one_pos.begin(), same_one_pos.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first < b.first)
                return true;
            else if (a.first > b.first)
                return false;
            else {
                return compSingBit(a, b);
            }
        });
        max_res_it = std::max_element(exist_zero_pos.begin(), exist_zero_pos.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first > b.first)
              return true;
            else if (a.first < b.first)
              return false;
            else {
              return compSingBit(a, b);
            }
        });
    } else if constexpr (std::is_same< decltype(bitwise_func), decltype(std::bit_or<bool>())>::value) {
        /*
         *  or:
         *  1 | 1 = 1, 1 | 0 = 1, 1 | x = 1
         *  0 | 0 = 0, 0 | x = x, x | x = x
         *   min: find '1' exist, the larger the better, others are '0'
         *   max: find the same '0', the larger the better, others are '1'
         * */
        min_res_it = std::max_element(exist_one_pos.begin(), exist_one_pos.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first < b.first)
                return true;
            else if (a.first > b.first)
                return false;
            else {
                return compSingBit(a, b);
            }
        });
        max_res_it = std::max_element(same_zero_pos.begin(), same_zero_pos.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first > b.first)
              return true;
            else if (a.first < b.first)
              return false;
            else {
              return compSingBit(a, b);
}
        });
    } else if constexpr (std::is_same< decltype(bitwise_func), decltype(std::bit_xor<bool>())>::value) {
        /*
         *  xor:
         *  0 ^ 0 = 0, 0 ^ 1 = 1, 1 ^ 1 = 0
         *  x ^ 0 = x, x ^ 1 = x, x ^ x = x
         *   min: find the different '0' + different '1', the larger the better, others are '0'
         *   max: find the same '0' + same '1', the larger the better, others are '1'
         * */
        min_res_it = std::max_element(union_min.begin(), union_min.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first < b.first)
              return true;
            else if (a.first > b.first)
              return false;
            else {
              return compSingBit(a, b);
            }
        });
        max_res_it = std::max_element(union_max.begin(), union_max.end(),
                                      [compSingBit](const signedBit& a, const signedBit& b) {
            if (a.first > b.first)
              return true;
            else if (a.first < b.first)
              return false;
            else {
              return compSingBit(a, b);
}
        });
    } else {
        assert(false && "unknown bit_wise operation");
    }

    std::string min_res_str(bit_size, '0'), max_res_str(bit_size, '1');
    min_res_str[0] = min_res_it->first ? '1' : '0';
    for (const auto& bit : min_res_it->second) {
        if (bit != SIZE_MAX) {
            min_res_str[bit] = '1';
        }
    }
    max_res_str[0] = max_res_it->first ? '1' : '0';
    for (const auto& bit : max_res_it->second) {
        if (bit != SIZE_MAX) {
            max_res_str[bit] = '0';
        }
    }

    /*
     * 5. convert binary string to decimal int64_t
     * */
    int64_t min_res = convertBin2Dec(min_res_str, bit_size);
    int64_t max_res = convertBin2Dec(max_res_str, bit_size);
    return std::make_pair(min_res, max_res);
}

void
checkBitWiseRes(const int64_t lhs_low, const int64_t lhs_high,
            const int64_t rhs_low, const int64_t rhs_high,
            std::ofstream& ofs)
{
    int64_t and_low = lhs_low & rhs_low;
    int64_t and_high = lhs_low & rhs_low;
    int64_t or_low = lhs_low | rhs_low;
    int64_t or_high = lhs_low | rhs_low;
    int64_t xor_low = lhs_low ^ rhs_low;
    int64_t xor_high = lhs_low ^ rhs_low;
#pragma omp parallel for
    for (int64_t i = lhs_low; i <= lhs_high; i++)
    {
#pragma omp parallel for
        for (int64_t j = rhs_low; j <= rhs_high; j++)
        {
            int64_t tmp_and = i & j;
            if (clearly)
                std::cout << "tmp_and: " << i << " & " << j << " = " << tmp_and << std::endl;
            and_low = tmp_and < and_low ? tmp_and : and_low;
            and_high = tmp_and > and_high ? tmp_and : and_high;

            int64_t tmp_or = i | j;
            if (clearly)
                std::cout << "tmp_or: " << i << " | " << j << " = " << tmp_or << std::endl;
            or_low = tmp_or < or_low ? tmp_or : or_low;
            or_high = tmp_or > or_high ? tmp_or : or_high;

            int64_t tmp_xor = i ^ j;
            if (clearly)
                std::cout << "tmp_xor: " << i << " ^ " << j << " = " << tmp_xor << std::endl;
            xor_low = tmp_xor < xor_low ? tmp_xor : xor_low;
            xor_high = tmp_xor > xor_high ? tmp_xor : xor_high;
        }
    }
    std::pair<int64_t, int64_t> and_res = bitwiseInterval(lhs_low, lhs_high, rhs_low, rhs_high,
                                                          std::bit_and<bool>());
    std::pair<int64_t, int64_t> or_res = bitwiseInterval(lhs_low, lhs_high, rhs_low, rhs_high,
                                                         std::bit_or<bool>());
    std::pair<int64_t, int64_t> xor_res = bitwiseInterval(lhs_low, lhs_high, rhs_low, rhs_high,
                                                         std::bit_xor<bool>());

    if (and_res.first == and_low && and_res.second == and_high) {
        printf("compute success with[%ld, %ld] and [%ld, %ld]. real range is [%ld, %ld]\n",
               lhs_low, lhs_high, rhs_low, rhs_high, and_low, and_high);
        ofs << "compute success with [" << lhs_low << ", " << lhs_high << "] and [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << and_low << ", " << and_high << "]" << std::endl;
    }
    else
    {
        printf("compute error with [%ld, %ld] and [%ld, %ld]. real range is [%ld, %ld], "
               "compute range is [%ld, %ld]\n", lhs_low, lhs_high, rhs_low, rhs_high, and_low,
               and_high, and_res.first, and_res.second);
        ofs << "compute error with [" << lhs_low << ", " << lhs_high << "] and [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << and_low << ", " << and_high << "], compute range is ["
            << and_res.first << ", " << and_res.second << "]" << std::endl;
    }

    if (or_res.first == or_low && or_res.second == or_high) {
        printf("compute success with[%ld, %ld] or [%ld, %ld]. real range is [%ld, %ld]\n",
               lhs_low, lhs_high, rhs_low, rhs_high, or_low, or_high);
        ofs << "compute success with [" << lhs_low << ", " << lhs_high << "] or [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << or_low << ", " << or_high << "]" << std::endl;
    }
    else
    {
        printf("compute error with [%ld, %ld] or [%ld, %ld]. real range is [%ld, %ld], "
               "compute range is [%ld, %ld]\n", lhs_low, lhs_high, rhs_low, rhs_high, or_low,
               or_high, or_res.first, or_res.second);
        ofs << "compute error with [" << lhs_low << ", " << lhs_high << "] or [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << or_low << ", " << or_high << "], compute range is ["
            << or_res.first << ", " << or_res.second << "]" << std::endl;
    }

    if (xor_res.first == xor_low && xor_res.second == xor_high) {
        printf("compute success with[%ld, %ld] xor [%ld, %ld]. real range is [%ld, %ld]\n",
               lhs_low, lhs_high, rhs_low, rhs_high, xor_low, xor_high);
        ofs << "compute success with [" << lhs_low << ", " << lhs_high << "] xor [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << xor_low << ", " << xor_high << "]" << std::endl;
    }
    else
    {
        printf("compute error with [%ld, %ld] xor [%ld, %ld]. real range is [%ld, %ld], "
               "compute range is [%ld, %ld]\n", lhs_low, lhs_high, rhs_low, rhs_high, xor_low,
               xor_high, xor_res.first, xor_res.second);
        ofs << "compute error with [" << lhs_low << ", " << lhs_high << "] xor [" << rhs_low << ", "
            << rhs_high << "]. real range is [" << xor_low << ", " << xor_high << "], compute range is ["
            << xor_res.first << ", " << xor_res.second << "]" << std::endl;
    }
}

int
main(int argc, char ** argv)
{
	char *	pEnd;
	int64_t lhs_low, lhs_high, rhs_low, rhs_high;
	bool	full_case_test = false;
	if (argc == 5)
	{
		lhs_low	 = strtol(argv[1], &pEnd, 10);
		lhs_high = strtol(argv[2], &pEnd, 10);
		rhs_low	 = strtol(argv[3], &pEnd, 10);
		rhs_high = strtol(argv[4], &pEnd, 10);
		assert(lhs_low <= lhs_high && rhs_low <= rhs_high);
		std::cout << "lhs low = " << lhs_low << ", lhs high = " << lhs_high
			  << ", rhs low = " << rhs_low << ", rhs high = " << rhs_high << std::endl;
	}
	else
	{
		std::cout << "full case test!" << std::endl;
		full_case_test = true;
	}

    std::ofstream ofs("mod_result");
    assert(ofs.is_open() && "error opening mod_result");

	if (!full_case_test)
	{
        clearly = true;
//		checkModRes(lhs_low, lhs_high, rhs_low, rhs_high, ofs);
        checkBitWiseRes(lhs_low, lhs_high, rhs_low, rhs_high, ofs);
	}
	else
	{
        const int64_t start_num = -5;
        const int64_t end_num = 5;
        #pragma omp parallel for
		for (int64_t fc_lhs_high = start_num; fc_lhs_high < end_num; fc_lhs_high++)
		{
            #pragma omp parallel for
			for (int64_t fc_lhs_low = start_num; fc_lhs_low < fc_lhs_high; fc_lhs_low++)
			{
                #pragma omp parallel for
				for (int64_t fc_rhs_high = start_num; fc_rhs_high < end_num; fc_rhs_high++)
				{
                    #pragma omp parallel for
					for (int64_t fc_rhs_low = start_num; fc_rhs_low < fc_rhs_high; fc_rhs_low++)
					{
//						checkModRes(fc_lhs_low, fc_lhs_high, fc_rhs_low, fc_rhs_high, ofs);
                        checkBitWiseRes(fc_lhs_low, fc_lhs_high, fc_rhs_low, fc_rhs_high, ofs);
					}
				}
			}
		}
	}
    ofs.close();
	return 0;
}
