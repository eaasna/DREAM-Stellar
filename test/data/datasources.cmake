cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/example_data.tar.gz
# declare_datasource (FILE example_data.tar.gz
#                     URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
#                     URL_HASH SHA256=7c2e7bdbf573cfe2314c8255080a5d57d966722e8bfb53712e03d87ce463ff15)


# declare_datasource (FILE expected_results.tar.gz
#                     URL https://ftp.imp.fu-berlin.de/pub/seiler/raptor/expected_results.tar.gz
#                     URL_HASH SHA256=2685ef95ebea074514f4736888b493857f0327514684ef88d798b3f25df5fd5a)

declare_datasource (FILE bin_0.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_0.fasta
                    URL_HASH SHA256=9a4a626225fdc66cfe2e05176adb9f140ea1d116d554463bc4163af5ad40a5b9)
declare_datasource (FILE bin_1.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_1.fasta
                    URL_HASH SHA256=3003ce921b4b232061964da79c54616c94e9f785628e0bac741510427c9d3760)
declare_datasource (FILE bin_2.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_2.fasta
                    URL_HASH SHA256=ef588e887403692bd770b2bc1f02eb66960552c963ea1ca1828356b8bbcbc38f)
declare_datasource (FILE bin_3.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_3.fasta
		            URL_HASH SHA256=e1934131d117dc9a12ebef54c180eab7b580648265ad56928d8570a33bf1daa1)
declare_datasource (FILE bin_4.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_4.fasta
                    URL_HASH SHA256=93f5a43a60c7beb175281f2f50674aa54e43320af8ac27976de8f96890c6d83b)
declare_datasource (FILE bin_5.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_5.fasta
                    URL_HASH SHA256=84baf896ea1f289384dde1319c512749674bcd0e5f61d92522f768f15189d3f3)
declare_datasource (FILE bin_6.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_6.fasta
                    URL_HASH SHA256=deb8c2cd17a59ec58d96c5ee9063b553917461661263f6ebeb44625a29076b87)
declare_datasource (FILE bin_7.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/bins/bin_7.fasta
                    URL_HASH SHA256=13fbd6105b6544df86ae948b0dc85aed2d6e8fb711942fb5cc283602c9ccbc26)

declare_datasource (FILE query.fq
                    URL ${CMAKE_SOURCE_DIR}/test/data/query.fq
                    URL_HASH SHA256=9eea085c8c36c6d1ff2d33ec6affe7382c78cc33c415072ad414c90060ed05d9)

declare_datasource (FILE various_chromosome_lengths.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/various_chromosome_lengths.fasta
                    URL_HASH SHA256=a6dccb2c5f434bea9a7afce246e2f50a6f59f96eeca126284274ea1581b7249a)
declare_datasource (FILE chromosome_metadata.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/chromosome_metadata.txt
                    URL_HASH SHA256=435e36bbdb0d481070aca7a86d3b264359bc0c5029fbdaa48f17b9e0f8ad306b)
declare_datasource (FILE 0overlap4bins.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/0overlap4bins.txt
                    URL_HASH SHA256=88ee7061dc9dd86e7bd4cd8e37424bd73c3a74a7574d0b93d5519580d7cf1470)
declare_datasource (FILE 20overlap4bins.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/20overlap4bins.txt
                    URL_HASH SHA256=88ee7061dc9dd86e7bd4cd8e37424bd73c3a74a7574d0b93d5519580d7cf1470)
declare_datasource (FILE 0overlap16bins.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/0overlap16bins.txt
                    URL_HASH SHA256=39b3c960ce05e23d3208ab766388819a1e82b8036ae303e0eae37b425f66ef36)
declare_datasource (FILE 20overlap16bins.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/split/20overlap16bins.txt
                    URL_HASH SHA256=b4d1c624381a18f120ecda7c4b7330fe34051503119ddefd8e08732fc20dca9f)

declare_datasource (FILE 8bins20window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window.ibf
                    URL_HASH SHA256=5ba1b39ff4a7c9e820b3cf738116d507f65c5d50c91221d77ff173d4297cb13b)
declare_datasource (FILE 8bins23window.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window.ibf
                    URL_HASH SHA256=2c528498ef692dde084c9d26b0583eb5bb74d640a1b777193764c3d667014c3d)

declare_datasource (FILE 8bins20window0error50pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window0error50pattern1overlap.out
                    URL_HASH SHA256=fb9e916a394acb4d6bcd53556cffe7300bee4915ed61d67a8abdc5d1f412d8a8)
declare_datasource (FILE 8bins20window0error50pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window0error50pattern40overlap.out
                    URL_HASH SHA256=880e7c2c50943436f96dc026f03d5aed89672df8a3b0c02b6b72d20105024f77)
declare_datasource (FILE 8bins20window0error100pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window0error100pattern1overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins20window0error100pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window0error100pattern40overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins20window1error50pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window1error50pattern1overlap.out
                    URL_HASH SHA256=55550efe40e0b9ba12c55be2f6ea6756279e249f2eb94a3474736377461b6b6f)
declare_datasource (FILE 8bins20window1error50pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window1error50pattern40overlap.out
                    URL_HASH SHA256=55550efe40e0b9ba12c55be2f6ea6756279e249f2eb94a3474736377461b6b6f)
declare_datasource (FILE 8bins20window1error100pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window1error100pattern1overlap.out
                    URL_HASH SHA256=29ad8d60ab1b3f12efd7e47a79f85dba36808e213f90c74236810b6801b32886)
declare_datasource (FILE 8bins20window1error100pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins20window1error100pattern40overlap.out
                    URL_HASH SHA256=29ad8d60ab1b3f12efd7e47a79f85dba36808e213f90c74236810b6801b32886)
declare_datasource (FILE 8bins23window0error50pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window0error50pattern1overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins23window0error50pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window0error50pattern40overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins23window0error100pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window0error100pattern1overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins23window0error100pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window0error100pattern40overlap.out
                    URL_HASH SHA256=6631c798ab13356c18938aca418f16e474226ad50e4e08bbf364cfc422ae0fe7)
declare_datasource (FILE 8bins23window1error50pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window1error50pattern1overlap.out
                    URL_HASH SHA256=55550efe40e0b9ba12c55be2f6ea6756279e249f2eb94a3474736377461b6b6f)
declare_datasource (FILE 8bins23window1error50pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window1error50pattern40overlap.out
                    URL_HASH SHA256=55550efe40e0b9ba12c55be2f6ea6756279e249f2eb94a3474736377461b6b6f)
declare_datasource (FILE 8bins23window1error100pattern1overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window1error100pattern1overlap.out
                    URL_HASH SHA256=f3e73e1b7d2826753d99f3ad1f48925e06ee8d744d00375c7fafaf9a123cdff1)
declare_datasource (FILE 8bins23window1error100pattern40overlap.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/8bins23window1error100pattern40overlap.out
                    URL_HASH SHA256=f3e73e1b7d2826753d99f3ad1f48925e06ee8d744d00375c7fafaf9a123cdff1)
