def calc_n(self, i):  # electrical coupling on
    i = g_func.toarray(i, self.stack.cell_numb, self.stack.cell.cathode.channel.nodes)
    self.stack.set_i(i)
    var0 = self.stack.cell.cathode.tar_cd * self.stack.cell.cathode.channel.nodes
    self.stack.i[int(self.stack.cell_numb / 2), -1] = var0 \
                                                      - sum(self.stack.i[int(self.stack.cell_numb / 2)]) \
                                                      + self.stack.i[int(self.stack.cell_numb / 2), -1]
    self.stack.update()
    var1 = np.matmul(self.b, self.stack.v)
    var2 = self.stack.resistance * np.matmul(self.c, self.stack.i.flatten(order='C'))
    return var1 - var2


def calc_solve(self, i):
    i = g_func.toarray(i, self.stack.cell_numb, self.stack.cell.cathode.channel.nodes)
    self.stack.set_i(i)
    self.stack.update()
    return self.stack.i.flatten(order='C') - self.stack.cell.cathode.tar_cd \
           * self.stack.cell_numb * self.stack.cell.cathode.channel.nodes


def calc_g(self, i):
    self.stack.set_i(g_func.toarray(i, self.stack.cell_numb, self.stack.cell.cathode.channel.nodes))
    var0 = self.stack.cell.cathode.tar_cd * self.stack.cell.cathode.channel.nodes
    self.stack.i[int(self.stack.cell_numb / 2), -1] = var0 \
                                                      - sum(self.stack.i[int(self.stack.cell_numb / 2)]) \
                                                      + self.stack.i[int(self.stack.cell_numb / 2), 0]
    self.stack.update()
    self.s = np.diag(self.stack.dv)
    return np.matmul(self.b, self.s) - self.stack.resistance * self.c


def calc_i(self):
    # fsolve(self.calc_n,self.stack.i.flatten(order='C'),
    # xtol= self.k_it,maxfev = self.max_it)
    # fsolve(self.calc_solve,self.stack.i.flatten(order='C'),
    # xtol= self.k_it,maxfev = self.max_it)
    # i =  fsolve(self.calc_n,self.stack.i.flatten(order='C'),
    #  fprime= self.calc_g,xtol= self.k_it,maxfev = self.max_it)
    # i = root(self.calc_n, self.stack.i.flatten(order='C'),
    # method='krylov', tol=self.k_it)
    self.stack.set_i(g_func.toarray(i.x,
                                    self.stack.cell_numb,
                                    self.stack.cell.cathode.channel.nodes))
    self.stack.update()
    self.calc_initial_current_density()