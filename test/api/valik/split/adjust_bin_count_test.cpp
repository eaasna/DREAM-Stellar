#include <gtest/gtest.h>

#include <valik/shared.hpp>

TEST(adjust_test, try_count)
{
    std::vector<size_t> bin_count_in{8, 63, 64, 65, 96, 97, 159, 2060, 4111};
    for (size_t b : bin_count_in)
    {
        size_t adjusted_count = valik::adjust_bin_count(b);

        if (b == std::numeric_limits<uint32_t>::max())
            EXPECT_EQ(adjusted_count, 64UL);    //default
        else if (b < 97UL)
            EXPECT_EQ(adjusted_count, 64UL);
        else if (b < 160UL)
            EXPECT_EQ(adjusted_count, 128UL);
        else
            EXPECT_TRUE(std::abs((int64_t)b - (int64_t)adjusted_count) <= 64);
    }
}
