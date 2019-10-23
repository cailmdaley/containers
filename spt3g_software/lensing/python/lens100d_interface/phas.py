import os
import numpy as np
import pickle as pk
from spt3g import core
from .. import utils


class RandomPhase(object):
    """
    Class for generating and caching a random state.
    """
    def __init__(
        self, size, lib_dir, seed=None, random=np.random.standard_normal
    ):
        self.size = size
        self.lib_dir = lib_dir
        self.random = random

        if True:
            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if not os.path.exists(lib_dir + "/state_%04d.pk" % 0):
                if seed != None:
                    np.random.seed(seed)
                else:
                    input(
                        "ds::sims::phas: press enter to initialize random seed."
                    )
                    np.random.seed()

                # store the internal state
                pk.dump(
                    np.random.get_state(),
                    open(lib_dir + "/state_%04d.pk" % 0, "wb"),
                )

            if not os.path.exists(lib_dir + "/sim_hash.pk"):
                pk.dump(self.hashdict(), open(lib_dir + "/sim_hash.pk", "wb"))

        utils.hash_check(
            pk.load(open(lib_dir + "/sim_hash.pk", "rb")), self.hashdict()
        )

        if seed != None:
            np.random.seed(seed)
            self.random(size=self.size)
            self.check_state_final(0)

    def hashdict(self):
        return {
            "size": self.size,
            "init": dict(
                list(
                    zip(
                        np.arange(0, 4),
                        pk.load(
                            open(self.lib_dir + "/state_%04d.pk" % 0, "rb"),
                            encoding="latin1",
                        ),
                    )
                )
            ),
        }

    def set_state(self, idx):
        np.random.set_state(self.get_state(idx))

    def get_state(self, idx):
        if not os.path.exists(self.lib_dir + "/state_%04d.pk" % idx):
            state = self.get_state(idx - 1)
            np.random.set_state(state)
            core.log_info("ds::sims::phas: caching state %d" % idx, unit="RandomPhase")
            self.random(size=self.size)
            pk.dump(
                np.random.get_state(),
                open(self.lib_dir + "/state_%04d.pk" % idx, "wb"),
            )
        return pk.load(
            open(self.lib_dir + "/state_%04d.pk" % idx, "rb"),
            encoding="latin1",
        )

    # Check the state after getting random numbers for set [idx]
    def check_state_final(self, idx):
        fstate = np.random.get_state()
        mstate = self.get_state(idx + 1)
        assert np.all([np.all(f == m) for f, m in zip(fstate, mstate)])
