#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "extmem.h"
/*线性搜索的关系选择算法*/
/*end_add:关系的结束块地址, next_addr:后继块的地址,point_addr:写入ַ的磁盘块号*/
int Linar_Search(int next_addr, Buffer *buf, int end_add, int point_addr)
{
    int i = 0, X = -1, Y = -1;
    int select_res = 0, p = 0;
    char str[5] = {0}, str1[5] = {0};
    unsigned char  *blk_io, *blk;
    blk_io = getNewBlockInBuffer(buf);
    while (next_addr <= end_add)
    {
        if (buf->numFreeBlk == 0)
        {
            freeBlockInBuffer(blk, buf);
            //freeBlockInBuffer(blk - ((next_addr - 2) * 64), &buf);
        }
        if ((blk = readBlockFromDisk(next_addr, buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        printf("读入磁盘块%d\n", next_addr);
        for (i = 0; i < 7; ++i)
        {
            for (int k = 0; k < 4; ++k)
            {
                str[k] = *(blk + i*8 + k);
            }
            X = atoi(str);
            for (int k = 0; k < 4; ++k)
            {
                str[k] = *(blk + i*8 + 4 + k);
            }
            Y = atoi(str);
            if (X == 30 && end_add ==16 || X == 23 && end_add == 48)
            {
                printf("(X=%d, Y=%d)\n", X, Y);
                p = select_res % 7;
                select_res += 1;
                itoa(X, str1, 10);
                for (int j = 0; j < 4; ++j)
                {
                    *(blk_io + 8*p + j) = str1[j];
                }
                itoa(Y, str1, 10);
                for (int j = 0; j < 4; ++j)
                {
                    *(blk_io + 8*p + j + 4) = str1[j];
                }
                if (p == 6)
                {
                    itoa(point_addr+1, str1, 10);
                    for (int j = 0; j < 4; ++j)
                    {
                        *(blk_io + 8*(p+1) + j) = str1[j];
                    }
                    writeBlockToDisk(blk_io, point_addr, buf);
                    blk_io = getNewBlockInBuffer(buf);
                    point_addr += 1;

                }
            }
        }
        for (int k = 0; k < 4; ++k)
        {
            str[k] = *(blk + i*8 + k);
        }
        next_addr = atoi(str);

    }
    if (p < 6)
    {
        itoa(point_addr+1, str1, 10);
        for (int i = 0; i < 4; ++i)
        {
            *(blk_io + 56 + i) = str1[i];
        }
        writeBlockToDisk(blk_io, point_addr, buf);
    }


    freeBlockInBuffer(blk_io, buf);
    printf("\n满足条件的元组一共%d个.\n", select_res);
    printf("IO读写一共%lu次.\n", buf->numIO);
    return next_addr;

}
/*冒泡排序,内排序*/
/*N是关系的总块数,m是子集合的个数,num代表第几个子集合,w_be是写入磁盘块的起始地址,r_be是读磁盘块的起始地址*/
void Bubble_Sort(int N, int m, Buffer *buf, int num, int w_be, int r_be)
{
    int X1 = -1, X2 = -2;
    unsigned char *blk, *blk0;
    char str1[5] = {0}, str2[5] = {0}, str3[5] = {0}, str4[5] = {0};

    char str[5] = {0};

    for (int j = 0; j < N / m; ++j)
    {
        if ((blk = readBlockFromDisk(j + N/m*(num-1) + r_be, buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        if (j == 0) blk0 = blk;
    }

    for (int k = 0; k < N / m * 8; ++k)
    {
        for (int l = 0; l < N / m * 8 - k - 1 - 1; ++l)
        {
            if (l == -1 + 8*(l/8+1)) continue;
            for (int o = 0; o < 4; ++o)
            {
                str1[o] = *(blk0 + 8*(l%8)+65*(l/8)+o);//X1
                str2[o] = *(blk0 + 8*(l%8)+65*(l/8)+o+4);//Y1
            }
            X1 = atoi(str1);
            for (int o = 0; o < 4; ++o)
            {
                if ((l + 1) == -1 + 8*((l+1)/8+1))
                {
                    str3[o] = *(blk0 + 8*((l+1+1)%8)+65*((l+1+1)/8)+o);//X2
                    str4[o] = *(blk0 + 8*((l+1+1)%8)+65*((l+1+1)/8)+o+4);//Y2
                }
                else
                {
                    str3[o] = *(blk0 + 8*((l+1)%8)+65*((l+1)/8)+o);//X2
                    str4[o] = *(blk0 + 8*((l+1)%8)+65*((l+1)/8)+o+4);//Y2
                }
            }
            X2 = atoi(str3);
            if (X1 > X2)
            {
                for (int o = 0; o < 4; ++o)
                {
                    if ((l + 1) == -1 + 8*((l+1)/8+1))
                    {
                        *(blk0 + 8*((l+1+1)%8)+65*((l+1+1)/8)+o) = str1[o];
                        *(blk0 + 8*((l+1+1)%8)+65*((l+1+1)/8)+o+4) = str2[o];
                        *(blk0 + 8*(l%8)+65*(l/8)+o) = str3[o];
                        *(blk0 + 8*(l%8)+65*(l/8)+o+4) = str4[o];
                    }
                    else
                    {
                        *(blk0 + 8*((l+1)%8)+65*((l+1)/8)+o) = str1[o];
                        *(blk0 + 8*((l+1)%8)+65*((l+1)/8)+o+4) = str2[o];
                        *(blk0 + 8*(l%8)+65*(l/8)+o) = str3[o];
                        *(blk0 + 8*(l%8)+65*(l/8)+o+4) = str4[o];
                    }

                }
            }
        }
    }
    for (int i = 0; i < N / m; ++i)
    {
        itoa(w_be+i + N/m*(num-1)+1, str, 10);
        for (int j = 0; j < 4; ++j)
            *(blk0 + 56 + j) = str[j];
        writeBlockToDisk(blk0+(i*65), w_be+i + N/m*(num-1), buf);
    }
    return;

}

/*TPMMS第二趟,归并排序*/
/*n是每个子集合的块数,r_be读取磁盘块的起始地址,w_be是写入磁盘块的起始地址,m子集合的个数*/
void Merge_Sort(int n, int r_be, int m, int w_be, Buffer *buf)
{
    unsigned char *blk;
    int X[5];
    int p[5] = {0};/*指向用了多少块*/
    int q[5] = {0};/*指向每个块用了多少元组*/
    int flag2[5] = {0};/*标记*/
    char str1[5] = {0};
    char str[5][5] = {0};
    unsigned char *blk_tmp[5] = {NULL};
    int num_in_blk = 0;
    int record = -1;

    memset(q, 0, sizeof(q));
    memset(flag2, 0, sizeof(flag2));
    num_in_blk = 0;

    for (int j = 0; j < m; ++j)
    {

        if ((blk = readBlockFromDisk(r_be+j*n+p[j], buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        p[j] += 1;
        blk_tmp[j] = blk;
    }
    blk = getNewBlockInBuffer(buf);

    while (p[0] <= n && p[1] <= n && p[2] <= n && p[3] <= n)
    {
        for (int i = 0; i < 4; ++i)
        {
            if (q[i] == 7 && p[i] < n)
            {
                freeBlockInBuffer(blk_tmp[i], buf);

                if ((blk_tmp[i] = readBlockFromDisk(r_be+p[i]+n*i, buf)) == NULL)
                {
                    perror("Reading Block Failed!\n");
                    return -1;
                }
                p[i] += 1;
                q[i] = 0;
            }

        }
        for (int i = 0; i < 4; ++i)
        {
            if (p[i] == n && q[i] == 7) flag2[i] = 1;
        }

        for (int j = 0; j < 4; ++j)
        {
            if (!flag2[0]) str[0][j] = *(blk_tmp[0] + j + 8*q[0]);//X1
            if (!flag2[1]) str[1][j] = *(blk_tmp[1] + j + 8*q[1]);//X2
            if (!flag2[2]) str[2][j] = *(blk_tmp[2] + j + 8*q[2]);//X3
            if (!flag2[3]) str[3][j] = *(blk_tmp[3] + j + 8*q[3]);//X4
        }
        for (int i = 0; i < 4; ++i)
        {
            if (!flag2[i]) X[i] = atoi(str[i]);
        }

        int Min = 1111111;
        int flag1 = 0;
        for (int k = 0; k < 4; ++k)
        {
            if (flag2[k]) continue;
            if(X[k] < Min)
            {
                Min = X[k];
                flag1 = k;
            }
        }
        for (int k = 0; k < 4; ++k)
        {
            *(blk + k + 8*num_in_blk) = str[flag1][k];
        }
        for (int k = 0; k < 4; ++k)
        {
            *(blk + k + 8*num_in_blk+4) = *(blk_tmp[flag1] + k + 4 + 8*q[flag1]);
        }
        num_in_blk += 1;
        if (num_in_blk == 7) {
            record += 1;
            num_in_blk = 0;
            itoa(w_be + record + 1, str1, 10);
            for (int y = 0; y < 4; ++y)
            {
                *(blk + 56 + y) = str1[y];
            }
            writeBlockToDisk(blk, w_be + record, buf);
            blk = getNewBlockInBuffer(buf);
        }
        q[flag1] += 1;

        if (p[0] == n && p[1] == n && p[2] == n && p[3] == n && q[0] == 7 && q[1] == 7 && q[2] == 7 && q[3] == 7) break;
    }

    return;
}

/*建立索引块*/
/*total_number:R/S的总块数,read_be:读磁盘的起始地址,write_be写磁盘的起始地址*/
void Create_Index(Buffer *buf, int total_number, int read_be, int write_be)
{
    unsigned char *blk, *blk0;
    int num = 0;
    int j = 0;
    char str[5] = {0}, str1[5] = {0};
    blk0 = getNewBlockInBuffer(buf);

    if (total_number % 7 == 0)
    {
        num = 1;
    }
    else
    {
        num = 0;
    }
    for (j = 0; j < total_number; ++j)
    {
        if ((blk = readBlockFromDisk(read_be+j, buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        for (int k = 0; k < 4; ++k)
        {
            str[k] = *(blk + k);
        }
        itoa(read_be+j,str1,10);
        for (int k = 0; k < 4; ++k)
        {
            *(blk0 + k + 8*(j%7)) = str[k];
            *(blk0 + k + 4 + 8*(j%7)) = str1[k];
        }
        if ((j+1) % 7 == 0)
        {
            memset(str1, 0, sizeof(str1));
            itoa(write_be + (j+1)/7, str1, 10);
            for (int k = 0; k < 4; ++k)
            {
                *(blk0 + 56 + k) = str1[k];
            }
            writeBlockToDisk(blk0, write_be-1+(j+1)/7, buf);
            freeBlockInBuffer(blk0, buf);
            blk0 = getNewBlockInBuffer(buf);

        }
        freeBlockInBuffer(blk, buf);

    }
    if (!num)
    {
        for (int i = total_number % 7; i < 8; ++i)
        {
            for (int j = 0; j < 8; ++j)
                *(blk0 + i*8+j) = '\0';
        }
        itoa(write_be+ j/7+1, str1, 10);
        for (int k = 0; k < 4 ;++k)
        {
            *(blk0 + 56 +k) = str1[k];
        }
        writeBlockToDisk(blk0, write_be+ j/7, buf);
    }
    return;
}

/*索引查找*/
/*total_number:R/S的总块数, read_be:读索引磁盘的起始地址ַ,write_be写磁盘的起始地址,key:键值*/
int Index_Search(int total_number, Buffer *buf, int read_be, int key, int write_be)
{
    int X[2] = {0}, Y[2] = {0}, X1 = -1, X2 = -2;
    char str[5] = {0}, str1[5] = {0}, str2[5] = {0}, str3[5] = {0};
    int record = 0,tmp = 0;
    unsigned char *blk1, *blk, *blkIO;
    int satisfy = 0;//满足条件的个数
    blkIO = getNewBlockInBuffer(buf);
    for (int i = 0; i < total_number / 7; ++i)
    {
        printf("读入索引块%d\n",read_be+i);
        if ((blk = readBlockFromDisk(read_be+i, buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        for (int j = 0; j < 4; ++j)
        {
            str2[j] = *(blk + 48 + j);
        }
        X2 = atoi(str2);//最大值
        for (int j = 0; j < 4; ++j)
        {
            str[j] = *(blk + j);
        }
        X1 = atoi(str);//最小值


        if (X1 > key || X2 < key)
        {
            printf("没有满足条件的元组.\n");
            freeBlockInBuffer(blk, buf);
        }
        else
        {
            for (int j = 0; j < 7; j++)
            {
                for (int k = 0; k < 4; ++k)
                {
                    str[k] = *(blk+k+8*j);
                    str1[k] = *(blk + k+8*j+4);
                }
                X[j%2] = atoi(str);
                Y[j%2] = atoi(str1);
                if (j >= 1)
                {
                    if (X[j%2] > key)
                    {
                        if (!record) {
                            j = 7;
                            break;
                        }
                        if (j % 2 == 0)//
                        {
                            if (X[j%2+1] <= key && j != 1) break;
                            printf("读入数据块%d\n", Y[j%2+1]);
                            if ((blk1 = readBlockFromDisk(Y[j%2+1], buf)) == NULL)
                            {
                                perror("Reading Block Failed!\n");
                                return -1;
                            }
                            for (int o = 0; o < 7; ++o)
                            {
                                for (int k = 0; k < 4; ++k)
                                {
                                    str[k] = *(blk1+k+8*o);
                                    str1[k] = *(blk1 + k+8*o+4);
                                }
                                X1 = atoi(str);
                                X2 = atoi(str1);
                                if (X1 == key)
                                {
                                    printf("(X=%d,Y=%d)\n", X1,X2);
                                    for (int o = 0; o < 4; ++o) {
                                        *(blkIO + o + (satisfy%7)*8) = str[o];
                                        *(blkIO + o + 4 + (satisfy%7)*8) = str1[o];
                                    }
                                    satisfy += 1;
                                    if (satisfy%7 == 0) {
                                        itoa(write_be+1, str3, 10);
                                        for (int o = 0; o < 4; ++o) {
                                            *(blkIO + 56 + o) = str3[o];
                                        }
                                        writeBlockToDisk(blkIO, write_be, buf);
                                        blkIO = getNewBlockInBuffer(buf);
                                        write_be += 1;
                                    }
                                }
                            }
                            freeBlockInBuffer(blk1, buf);
                        }
                        else
                        {
                            if (X[j%2-1] <= key && j != 1) break;
                            printf("读入数据块%d\n", Y[j%2-1]);
                            if ((blk1 = readBlockFromDisk(Y[j%2-1], buf)) == NULL)
                            {
                                perror("Reading Block Failed!\n");
                                return -1;
                            }
                            for (int o = 0; o < 7; ++o)
                            {
                                for (int k = 0; k < 4; ++k)
                                {
                                    str[k] = *(blk1+k+8*o);
                                    str1[k] = *(blk1 + k+8*o+4);
                                }
                                X1 = atoi(str);
                                X2 = atoi(str1);
                                if (X1 == key)
                                {
                                    printf("(X=%d,Y=%d)\n", X1,X2);
                                    for (int o = 0; o < 4; ++o) {
                                        *(blkIO + o + (satisfy%7)*8) = str[o];
                                        *(blkIO + o + 4 + (satisfy%7)*8) = str1[o];
                                    }
                                    satisfy += 1;
                                    if (satisfy%7 == 0) {
                                        itoa(write_be+1, str3, 10);
                                        for (int o = 0; o < 4; ++o) {
                                            *(blkIO + 56 + o) = str3[o];
                                        }
                                        writeBlockToDisk(blkIO, write_be, buf);
                                        blkIO = getNewBlockInBuffer(buf);
                                        write_be += 1;
                                    }
                                }

                            }
                            freeBlockInBuffer(blk1, buf);
                        }
                    }
                    else if (X[j%2] == key)
                    {
                        record += 1;
                        if (record == 1)//第一次进入
                        {
                            if (j & 1 == 0) {
                                tmp = j % 2 - 1;
                            } else {
                                tmp = j % 2 + 1;
                            }
                            for (int y = 0; y < 2; ++y)
                            {
                                printf("读入数据块%d\n", Y[(tmp+y)%2]);
                                if ((blk1 = readBlockFromDisk(Y[(tmp+y)%2], buf)) == NULL)
                                {
                                    perror("Reading Block Failed!\n");
                                    return -1;
                                }
                                for (int o = 0; o < 7; ++o)
                                {
                                    for (int k = 0; k < 4; ++k)
                                    {
                                        str[k] = *(blk1+k+8*o);
                                        str1[k] = *(blk1 + k+8*o+4);
                                    }
                                    X1 = atoi(str);
                                    X2 = atoi(str1);
                                    if (X1 == key)
                                    {
                                        printf("(X=%d,Y=%d)\n", X1,X2);
                                        for (int o = 0; o < 4; ++o) {
                                            *(blkIO + o + (satisfy%7)*8) = str[o];
                                            *(blkIO + o + 4 + (satisfy%7)*8) = str1[o];
                                        }
                                        satisfy += 1;
                                        if (satisfy%7 == 0) {
                                            itoa(write_be+1, str3, 10);
                                            for (int o = 0; o < 4; ++o) {
                                                *(blkIO + 56 + o) = str3[o];
                                            }
                                            writeBlockToDisk(blkIO, write_be, buf);
                                            blkIO = getNewBlockInBuffer(buf);
                                            write_be += 1;
                                        }
                                    }
                                }
                                freeBlockInBuffer(blk1, buf);
                            }
                        }
                        else
                        {
                            printf("读入数据块%d\n", Y[j%2]);
                            if ((blk1 = readBlockFromDisk(Y[j%2], buf)) == NULL)
                            {
                                perror("Reading Block Failed!\n");
                                return -1;
                            }
                            for (int o = 0; o < 7; ++o)
                            {
                                for (int k = 0; k < 4; ++k)
                                {
                                    str[k] = *(blk1+k+8*o);
                                    str1[k] = *(blk1 + k+8*o+4);
                                }
                                X1 = atoi(str);
                                X2 = atoi(str1);
                                if (X1 == key)
                                {
                                    printf("(X=%d,Y=%d)\n", X1,X2);
                                    for (int o = 0; o < 4; ++o) {
                                        *(blkIO + o + (satisfy%7)*8) = str[o];
                                        *(blkIO + o + 4 + (satisfy%7)*8) = str1[o];
                                    }
                                    satisfy += 1;
                                    if (satisfy%7 == 0) {
                                        itoa(write_be+1, str3, 10);
                                        for (int o = 0; o < 4; ++o) {
                                            *(blkIO + 56 + o) = str3[o];
                                        }
                                        writeBlockToDisk(blkIO, write_be, buf);
                                        blkIO = getNewBlockInBuffer(buf);
                                        write_be += 1;
                                    }
                                }

                            }
                            freeBlockInBuffer(blk1, buf);
                        }
                    }
                }
            }
        }
    }
    if (satisfy % 7 != 0) {
        itoa(write_be+1, str3, 10);
        for (int o = 0; o < 4; ++o) {
            *(blkIO + 56 + o) = str3[o];
        }
        writeBlockToDisk(blkIO, write_be, buf);
    }
    return satisfy;
}

/*实现关系投影并去重*/
/*end_add:关系的结束块地址ַ, next_addr:关系的后继地址ַ,point_addr:写入磁盘的结果地址*/
int Project_Based(int next_addr, Buffer *buf, int end_add, int point_addr)
{
    int i = 0, X = -1, Y = -1;
    int project_res = 0, p = 0;
    char str[5] = {0}, str1[5] = {0};
    int flag[70] = {0};
    unsigned char  *blk_io, *blk;
    blk_io = getNewBlockInBuffer(buf);
    while (next_addr <= end_add)
    {
        if (buf->numFreeBlk == 0)
        {
            freeBlockInBuffer(blk, buf);
            //freeBlockInBuffer(blk - ((next_addr - 2) * 64), &buf);
        }
        if ((blk = readBlockFromDisk(next_addr, buf)) == NULL)
        {
            perror("Reading Block Failed!\n");
            return -1;
        }
        printf("读入数据块%d\n", next_addr);
        for (i = 0; i < 7; ++i)
        {
            for (int k = 0; k < 4; ++k)
            {
                str[k] = *(blk + i*8 + k);
            }
            X = atoi(str);
            if (!flag[X]) {
                printf("(X=%d)\n", X);
                flag[X] = 1;
                p = project_res % 14;
                project_res += 1;
                itoa(X, str1, 10);
                for (int j = 0; j < 4; ++j)
                {
                    *(blk_io + 4*p + j) = str1[j];
                }

                if (p == 13)
                {
                    itoa(point_addr+1, str1, 10);
                    for (int j = 0; j < 4; ++j)
                    {
                        *(blk_io + 4*(p+1) + j) = str1[j];
                    }
                    writeBlockToDisk(blk_io, point_addr, buf);
                    blk_io = getNewBlockInBuffer(buf);
                    point_addr += 1;

                }
            }
        }
        for (int k = 0; k < 4; ++k)
        {
            str[k] = *(blk + i*8 + k);
        }
        next_addr = atoi(str);
    }
    if (p < 13)
    {
        itoa(point_addr+1, str1, 10);
        for (int i = 0; i < 4; ++i)
        {
            *(blk_io + 56 + i) = str1[i];
        }
        writeBlockToDisk(blk_io, point_addr, buf);
    }

    freeBlockInBuffer(blk_io, buf);

    return project_res;


}

/*基于排序的连接算法*/
/*R_readbe:R排序后的起始地址,S_readbe:S排序后的起始地址*/
/*n_R:R磁盘块的总块数.n_S:S磁盘块的总块数,write_be:写入的结果磁盘块的起始地址*/
int Relation_Connect(int R_readbe, int S_readbe, int n_R, int n_S, int write_be, Buffer *buf)
{
    int res = 0;
    unsigned char *blk;
    int X[5];
    int p[5] = {0};
    int q[5] = {0};
    int flag2[5] = {0};
    char str1[5] = {0};
    char str[5][5] = {0};
    char strb[5] = {0}, stre[5] = {0};
    unsigned char *blk_tmp[5] = {NULL};
    unsigned char *blk_t, *blk_t0;
    int num_in_blk = 0;
    int record = -1;
    int not_move1 = 0, not_move2 = 0; /*标记R的指针,标记S的指针ֵ*/
    int yy = -1, ff = 0, f_flag = 0;
    memset(q, 0, sizeof(q));
    memset(flag2, 0, sizeof(flag2));
    num_in_blk = 0;

    if ((blk_tmp[0] = readBlockFromDisk(R_readbe+p[0], buf)) == NULL)
    {
        perror("Reading Block Failed!\n");
        return -1;
    }
    p[0] += 1;
    if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1], buf)) == NULL)
    {
        perror("Reading Block Failed!\n");
        return -1;
    }
    p[1] += 1;
    blk_t0 = blk_tmp[1];
    blk_t = blk_tmp[1];
    blk = getNewBlockInBuffer(buf);

    while (p[0] <= n_R && p[1] <= n_S)
    {
        if (q[0] == 7 && p[0] < n_R)
        {
            freeBlockInBuffer(blk_tmp[0], buf);

            if ((blk_tmp[0] = readBlockFromDisk(R_readbe+p[0], buf)) == NULL)
            {
                perror("Reading Block Failed!\n");
                return -1;
            }
            p[0] += 1;
            q[0] = 0;
        }
        if (q[1] == 7 && p[1] < n_S)
        {
            freeBlockInBuffer(blk_tmp[1], buf);

            if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1], buf)) == NULL)
            {
                perror("Reading Block Failed!\n");
                return -1;
            }
            blk_t = blk_tmp[1];
            p[1] += 1;
            q[1] = 0;
        }
        if ((p[0] == n_R && q[0] == 7) || (p[0] == n_R && q[0] == 7)) break;
        blk_tmp[1] = blk_t;
        for (int j = 0; j < 4; ++j)
        {

            str[0][j] = *(blk_tmp[0] + j + 8*q[0]);
            str[2][j] = *(blk_tmp[0] + j + 8*q[0] + 4);
            str[1][j] = *(blk_tmp[1] + j + 8*q[1]);
            str[3][j] = *(blk_tmp[1] + j + 8*q[1] + 4);
        }
        for (int i = 0; i < 4; ++i)
        {
            X[i] = atoi(str[i]);
        }
        while(q[0] < 7 && q[1] < 7)
        {
            blk_tmp[1] = blk_t;

            if (X[0] < X[1]) {
                q[0] += 1;
                for (int o = 0; o < 4; ++o) {
                    str[0][o] = *(blk_tmp[0] + o + 8*q[0]);
                    str[2][o] = *(blk_tmp[0] + o + 8*q[0] + 4);
                }
                X[0] = atoi(str[0]);
            }
            else if (X[0] > X[1]) {

                if (!ff) {
                    q[1] += 1;
                } else {
                    ff = 0;
                }
                if (q[1] == 7) break;
                for (int o = 0; o < 4; ++o) {
                    str[1][o] = *(blk_tmp[1] + o + 8*q[1]);
                    str[3][o] = *(blk_tmp[1] + o + 8*q[1] + 4);
                }
                X[1] = atoi(str[1]);
            }
            else {
                blk_tmp[1] = blk_t;
                if (write_be == 710) {
                    printf("11\n");
                }
                ff = 1;
                strcpy(strb, str[1]);
                strcpy(stre, str[3]);
                not_move1 = q[0];
                not_move2 = q[1];
                while (X[0] == X[1] && not_move2 < 7) {
                    res += 1;
                    if (num_in_blk == 6) {

                        for (int o = 0; o < 4; ++o) {
                            *(blk + o + 8*num_in_blk) = str[0][o];
                            *(blk + o + 8*num_in_blk + 4) = str[2][o];
                        }
                        itoa(write_be+1, str1, 10);
                        for (int o = 0; o < 4; ++o) {
                            *(blk + o +56) = str1[o];
                        }
                        printf("注:结果写入磁盘:%d.\n", write_be);
                        writeBlockToDisk(blk, write_be, buf);
                        write_be += 1;
                        num_in_blk = 0;
                        blk = getNewBlockInBuffer(buf);
                        for (int o = 0; o < 4; ++o) {
                            *(blk + o + 8*num_in_blk) = str[1][o];
                            *(blk + o + 8*num_in_blk + 4) = str[3][o];
                        }
                        num_in_blk += 1;

                    }
                    else if (num_in_blk + 2<= 7) {
                        for (int o = 0; o < 4; ++o) {
                            *(blk + o + 8*num_in_blk) = str[0][o];
                            *(blk + o + 8*num_in_blk + 4) = str[2][o];
                            *(blk + o + 8*(num_in_blk+1)) = str[1][o];
                            *(blk + o + 8*(num_in_blk+1) + 4) = str[3][o];

                        }
                        num_in_blk += 2;
                        if (num_in_blk == 7) {
                            itoa(write_be+1, str1, 10);
                            for (int o = 0; o < 4; ++o) {
                                *(blk + o +56) = str1[o];
                            }
                            printf("注:结果写入磁盘:%d.\n", write_be);
                            writeBlockToDisk(blk, write_be, buf);
                            write_be += 1;
                            num_in_blk = 0;
                            blk = getNewBlockInBuffer(buf);
                        }
                    }

                    not_move2 += 1;
                    if (not_move2 < 7) {
                        for (int o = 0; o < 4; ++o) {
                            str[1][o] = *(blk_tmp[1] + o + 8*not_move2);
                            str[3][o] = *(blk_tmp[1] + o + 8*not_move2 + 4);
                        }
                        X[1] = atoi(str[1]);
                    }
                    if (not_move2 % 7 == 0) {
                        yy += 1;
                        f_flag += 1;
                        if (buf->numFreeBlk == 0 && f_flag != 1) {
                            freeBlockInBuffer(blk_tmp[1], buf);
                        }
                        if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1]+yy, buf)) == NULL)
                        {
                            perror("Reading Block Failed!\n");
                            return -1;
                        }
                        not_move2 = 0;
                        for (int o = 0; o < 4; ++o) {
                            str[1][o] = *(blk_tmp[1] + o + 8*not_move2);
                            str[3][o] = *(blk_tmp[1] + o + 8*not_move2 + 4);
                        }
                        X[1] = atoi(str[1]);
                    }

                }
                if (f_flag != 0) freeBlockInBuffer(blk_tmp[1], buf);
                X[1] = X[0];
                strcpy(str[1], strb);
                strcpy(str[3], stre);
                yy = -1;
                f_flag = 0;
                q[0] += 1;
                if (q[0] < 7) {
                    for (int o = 0; o < 4; ++o) {
                        str[0][o] = *(blk_tmp[0] + o + 8*q[0]);
                        str[2][o] = *(blk_tmp[0] + o + 8*q[0] + 4);
                    }
                    X[0] = atoi(str[0]);
                }
            }
        }
    }
    return res;

}

/*关系的交*/
/*R_readbe:R排序后的起始地址,S_readbe:S排序后的起始地址ַ*/
/*n_R:R磁盘块的总块数,n_S:S磁盘块的总块数,write_be:写入的结果磁盘块的起始地址*/
int Intersection(int R_readbe, int S_readbe, int n_R, int n_S, int write_be, Buffer *buf)
{
    int res = 0;
    unsigned char *blk;
    int X[5], Y[5];
    int p[5] = {0};
    int q[5] = {0};
    int flag2[5] = {0};
    char str1[5] = {0};
    char str[5][5] = {0};
    char strb[5] = {0}, stre[5] = {0};
    unsigned char *blk_tmp[5] = {NULL};
    unsigned char *blk_t, *blk_t0;
    int num_in_blk = 0;
    int record = -1;
    int not_move1 = 0, not_move2 = 0;
    int yy = -1, ff = 0, f_flag = 0;
    memset(q, 0, sizeof(q));
    memset(flag2, 0, sizeof(flag2));
    num_in_blk = 0;

    if ((blk_tmp[0] = readBlockFromDisk(R_readbe+p[0], buf)) == NULL)
    {
        perror("Reading Block Failed!\n");
        return -1;
    }
    p[0] += 1;

    if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1], buf)) == NULL)
    {
        perror("Reading Block Failed!\n");
        return -1;
    }
    p[1] += 1;
    blk_t0 = blk_tmp[1];
    blk_t = blk_tmp[1];
    blk = getNewBlockInBuffer(buf);

    while (p[0] <= n_R && p[1] <= n_S)
    {
        if (q[0] == 7 && p[0] < n_R)
        {
            freeBlockInBuffer(blk_tmp[0], buf);

            if ((blk_tmp[0] = readBlockFromDisk(R_readbe+p[0], buf)) == NULL)
            {
                perror("Reading Block Failed!\n");
                return -1;
            }
            p[0] += 1;
            q[0] = 0;
        }
        if (q[1] == 7 && p[1] < n_S)
        {
            freeBlockInBuffer(blk_tmp[1], buf);
            if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1], buf)) == NULL)
            {
                perror("Reading Block Failed!\n");
                return -1;
            }
            blk_t = blk_tmp[1];
            p[1] += 1;
            q[1] = 0;
        }
        if ((p[0] == n_R && q[0] == 7) || (p[0] == n_R && q[0] == 7)) break;
        blk_tmp[1] = blk_t;
        for (int j = 0; j < 4; ++j)
        {
            str[0][j] = *(blk_tmp[0] + j + 8*q[0]);
            str[2][j] = *(blk_tmp[0] + j + 8*q[0] + 4);
            str[1][j] = *(blk_tmp[1] + j + 8*q[1]);
            str[3][j] = *(blk_tmp[1] + j + 8*q[1] + 4);
        }
        X[0] = atoi(str[0]);
        Y[0] = atoi(str[2]);
        X[1] = atoi(str[1]);
        Y[1] = atoi(str[3]);
        while(q[0] < 7 && q[1] < 7)
        {
            blk_tmp[1] = blk_t;
            if (X[0] < X[1]) {
                q[0] += 1;
                for (int o = 0; o < 4; ++o) {
                    str[0][o] = *(blk_tmp[0] + o + 8*q[0]);
                    str[2][o] = *(blk_tmp[0] + o + 8*q[0] + 4);
                }
                X[0] = atoi(str[0]);
                Y[0] = atoi(str[2]);
            }
            else if (X[0] > X[1]) {
                if (!ff) {
                    q[1] += 1;
                } else {
                    ff = 0;
                }
                if (q[1] == 7) break;
                for (int o = 0; o < 4; ++o) {
                    str[1][o] = *(blk_tmp[1] + o + 8*q[1]);
                    str[3][o] = *(blk_tmp[1] + o + 8*q[1] + 4);
                }
                X[1] = atoi(str[1]);
                Y[1] = atoi(str[3]);
            }
            else {
                blk_tmp[1] = blk_t;
                if (write_be == 710) {
                    printf("11\n");
                }
                ff = 1;
                strcpy(strb, str[1]);
                strcpy(stre, str[3]);
                not_move1 = q[0];
                not_move2 = q[1];
                while (X[0] == X[1] && not_move2 < 7) {
                    if (X[0] == X[1] && Y[0] == Y[1]) {
                        printf("(X=%d,Y=%d)\n", X[0],Y[0]);
                        res += 1;
                        if (num_in_blk == 6) {
                            for (int o = 0; o < 4; ++o) {
                                *(blk + o + 8*num_in_blk) = str[0][o];
                                *(blk + o + 8*num_in_blk + 4) = str[2][o];
                            }
                            itoa(write_be+1, str1, 10);
                            for (int o = 0; o < 4; ++o) {
                                *(blk + o +56) = str1[o];
                            }
                            printf("注:结果写入磁盘:%d.\n", write_be);
                            writeBlockToDisk(blk, write_be, buf);
                            write_be += 1;
                            num_in_blk = 0;
                            blk = getNewBlockInBuffer(buf);
                        }
                        else if (num_in_blk < 6) {
                            for (int o = 0; o < 4; ++o) {
                                *(blk + o + 8*num_in_blk) = str[0][o];
                                *(blk + o + 8*num_in_blk + 4) = str[2][o];
                            }
                            num_in_blk += 1;
                        }
                        break;
                    }
                    not_move2 += 1;
                    if (not_move2 < 7) {
                        for (int o = 0; o < 4; ++o) {
                            str[1][o] = *(blk_tmp[1] + o + 8*not_move2);
                            str[3][o] = *(blk_tmp[1] + o + 8*not_move2 + 4);
                        }
                        X[1] = atoi(str[1]);
                        Y[1] = atoi(str[3]);
                    }
                    if (not_move2 % 7 == 0) {
                        yy += 1;
                        f_flag += 1;
                        if (buf->numFreeBlk == 0 && f_flag != 1) {
                            freeBlockInBuffer(blk_tmp[1], buf);
                        }
                        if ((blk_tmp[1] = readBlockFromDisk(S_readbe+p[1]+yy, buf)) == NULL)
                        {
                            perror("Reading Block Failed!\n");
                            return -1;
                        }
                        not_move2 = 0;
                        for (int o = 0; o < 4; ++o) {
                            str[1][o] = *(blk_tmp[1] + o + 8*not_move2);
                            str[3][o] = *(blk_tmp[1] + o + 8*not_move2 + 4);
                        }
                        X[1] = atoi(str[1]);
                        Y[1] = atoi(str[3]);
                    }

                }
                if (f_flag != 0) freeBlockInBuffer(blk_tmp[1], buf);
                X[1] = X[0];
                strcpy(str[1], strb);
                strcpy(str[3], stre);
                X[1] = atoi(str[1]);
                Y[1] = atoi(str[3]);
                yy = -1;
                f_flag = 0;
                q[0] += 1;
                if (q[0] < 7) {
                    for (int o = 0; o < 4; ++o) {
                        str[0][o] = *(blk_tmp[0] + o + 8*q[0]);
                        str[2][o] = *(blk_tmp[0] + o + 8*q[0] + 4);
                    }
                    X[0] = atoi(str[0]);
                    Y[0] = atoi(str[2]);
                }
            }
        }
    }
    if (num_in_blk > 0) {
        for (int o = 0; o < 4; ++o) {
            *(blk + o + 8*num_in_blk) = str[0][o];
            *(blk + o + 8*num_in_blk + 4) = str[2][o];
        }
        itoa(write_be+1, str1, 10);
        for (int o = 0; o < 4; ++o) {
            *(blk + o +56) = str1[o];
        }
        printf("注:结果写入磁盘:%d.\n", write_be);
        writeBlockToDisk(blk, write_be, buf);
        write_be += 1;
        num_in_blk = 0;
        blk = getNewBlockInBuffer(buf);
    }
    return res;
}

int main(int argc, char **argv)
{
    Buffer buf;
    unsigned char *blk, *blk0;
    char str[5] = {0}, str1[5] = {0}, str2[5] = {0};
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    /*线性搜索*/
    int next_addr = 1;
    /*R*/
    printf("基于线性搜索的选择算法R.A=30:\n\n");
    next_addr = Linar_Search(next_addr, &buf, 16, 100);
    printf("注:结果写入磁盘:100\n");
    freeBuffer(&buf);
//    if (!initBuffer(520, 64, &buf))
//    {
//        perror("Buffer Initialization Failed!\n");
//        return -1;
//    }
//    /*S*/
//    printf("\n基于线性搜索的选择算法S.C=23:\n\n");
//    Linar_Search(next_addr, &buf, 48, 103);
//    printf("注:结果写入磁盘:103\n");
//    freeBuffer(&buf);

    /*TPMMS*/
    /*R*/
    int m = 4;/*子集合个数*/
    int R = 16, S = 32;

    int X1 = -1, X2 = -1;
    /*R的第一趟*/
    for (int i = 1; i <= m; ++i)
    {
        if (!initBuffer(520, 64, &buf))
        {
            perror("Buffer Initialization Failed!\n");
            return -1;
        }
        Bubble_Sort(R, m, &buf, i, 110, 1);
        free(&buf);
    }
    /*R的第二趟*/
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    Merge_Sort(R/m,110,m,300, &buf);
    free(&buf);
    /*S的第一趟*/
    for (int i = 1; i <= m; ++i)
    {
        if (!initBuffer(520, 64, &buf))
        {
            perror("Buffer Initialization Failed!\n");
            return -1;
        }
        Bubble_Sort(S, m, &buf, i, 130, 17);
        free(&buf);
    }
    /*S的第二趟*/
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    Merge_Sort(S/m,130,m,400, &buf);
    free(&buf);
    printf("关系R排序后输出到文件300.blk到315.blk.\n");
    printf("关系S排序后输出到文件400.blk到431.blk.\n");

    /*索引的选择算法*/
    /*R*/
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }

    printf("\n基于索引的选择算法R.A=30:\n\n");

    /*1.350-352为R的索引块*/

    Create_Index(&buf, R, 300, 350);
    free(&buf);
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    /*2.R的索引查找结果放在:500*/
    int res = Index_Search(R, &buf, 350, 30, 500);
    printf("注:结果写入磁盘:500\n");
    printf("满足选择条件的元组一共%d个.\n", res);
    printf("IO读写一共%d次.\n",buf.numIO);
    free(&buf);
    /*实现关系投影算法*/
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    int res1 = Project_Based(300, &buf, 315, 600);
    printf("注:结果写入磁盘:600\n\n");

    printf("关系R上的A属性满足投影(去重)的属性值一共%d个.\n", res1);
    free(&buf);

    /*实现基于排序的连接算法*/
    if (!initBuffer(520, 64, &buf))_
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    int res2 = Relation_Connect(300, 400, R, S, 700, &buf);
    printf("总共连接%d次.\n", res2);
    free(&buf);
    /*实现关系的交运算*/
    if (!initBuffer(520, 64, &buf))
    {
        perror("Buffer Initialization Failed!\n");
        return -1;
    }
    int res3 = Intersection(300,400,R,S,800, &buf);
    printf("S和R的交集有%d个元组.\n", res3);
    free(&buf);

    return 0;

}
