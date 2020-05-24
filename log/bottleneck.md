#  Possible places of bottleneck according to perf
LLcol llcol = colspace[joffset];
Tval w = llcol.cval * colScale;
Tind j = llcol.row;
Tval f = w/wdeg;
Tind jhead = a.cols[j];
Tind khead = a.cols[k];

### LLcol llcol = colspace[joffset];
3.11 :	  429d0e:       mov    -0x44(%rbp),%eax
0.13 :	  429d11:       movslq %eax,%rdx
0.06 :	  429d14:       lea    -0x150(%rbp),%rax
0.01 :	  429d1b:       mov    %rdx,%rsi
0.50 :	  429d1e:       mov    %rax,%rdi
0.00 :	  429d21:       callq  42d0b6 <std::vector<LLcol, std::allocator<LLcol>>::operator[](unsigned long)>
0.00 :	  429d26:       mov    0x8(%rax),%rdx
3.03 :	  429d2a:       mov    (%rax),%rax
0.30 :	  429d2d:       mov    %rax,-0x1c0(%rbp)
0.57 :	  429d34:       mov    %rdx,-0x1b8(%rbp)

### Tval w = llcol.cval * colScale;
0.46 :	  429d3b:       vmovsd -0x1b8(%rbp),%xmm0
2.34 :	  429d43:       vmulsd -0x40(%rbp),%xmm0,%xmm0
2.36 :	  429d48:       vmovsd %xmm0,-0x60(%rbp)

### Tind j = llcol.row;
0.48 :	  429d4d:       mov    -0x1c0(%rbp),%eax
0.03 :	  429d53:       mov    %eax,-0x1c4(%rbp)

### Tval f = w/wdeg;
0.00 :	  429d59:       vmovsd -0x60(%rbp),%xmm0
2.13 :	  429d5e:       vdivsd -0x38(%rbp),%xmm0,%xmm0
9.47 :	  429d63:       vmovsd %xmm0,-0x1d0(%rbp)

### Tind jhead = a.cols[j];
2.42 :	  429ebb:       mov    -0x1c4(%rbp),%eax
0.97 :	  429ec1:       cltq
0.15 :	  429ec3:       mov    -0x1f0(%rbp),%rdx
0.02 :	  429eca:       add    $0x8,%rdx
0.29 :	  429ece:       mov    %rax,%rsi
0.00 :	  429ed1:       mov    %rdx,%rdi
0.06 :	  429ed4:       callq  41fdba <std::vector<int, std::allocator<int> >::operator[](unsigned long)>
0.18 :	  429ed9:       mov    (%rax),%eax
7.26 :	  429edb:       mov    %eax,-0x74(%rbp)

### Tind khead = a.cols[k];
2.97 :	  429f79:       mov    -0x68(%rbp),%eax
1.22 :	  429f7c:       cltq
0.25 :	  429f7e:       mov    -0x1f0(%rbp),%rdx
0.10 :	  429f85:       add    $0x8,%rdx
0.43 :	  429f89:       mov    %rax,%rsi
0.00 :	  429f8c:       mov    %rdx,%rdi
0.00 :	  429f8f:       callq  41fdba <std::vector<int, std::allocator<int> >::operator[](unsigned long)>
0.00 :	  429f94:       mov    (%rax),%eax
4.29 :	  429f96:       mov    %eax,-0x78(%rbp)



