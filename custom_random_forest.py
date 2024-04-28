from typing import Tuple, List, Any
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator
from multiprocessing import Pool


class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom implementation of a Random Forest classifier.
    """

    def __init__(
        self,
        n_estimators: int = 10,
        max_depth: int = None,
        max_features: int = None,
        random_state: int = 42,
    ):
        """
        Initializes the RandomForestClassifierCustom instance.

        Args:
            n_estimators (int): The number of trees in the forest.
            max_depth (int): The maximum depth of the trees.
            max_features (int): The maximum number of features considered for splitting a node.
            random_state (int): A seed used by the random number generator.
        """
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []  # List to store each tree in the forest.
        self.feat_ids_by_tree = []  # List to store feature indices for each tree.

    def _fit_tree(
        self, args: Tuple[Any, np.ndarray, int, int, int]
    ) -> Tuple[DecisionTreeClassifier, np.ndarray]:
        """
        Fits a single decision tree using a subset of data and features.

        Args:
            args (tuple): Contains the dataset, labels, and parameters for tree fitting.

        Returns:
            A tuple containing the fitted tree and the feature indices used.
        """
        X, y, max_features, max_depth, random_state = args
        np.random.seed(random_state)
        feat_ids = np.random.choice(range(X.shape[1]), size=max_features, replace=False)
        pseudo_ids = np.random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
        pseudo_X = X[pseudo_ids, :][:, feat_ids]
        pseudo_y = y[pseudo_ids]
        dt_clf = DecisionTreeClassifier(
            max_depth=max_depth, max_features=max_features, random_state=random_state
        )
        dt_clf.fit(pseudo_X, pseudo_y)
        return dt_clf, feat_ids

    def fit(
        self, X: np.ndarray, y: np.ndarray, n_jobs: int = 1
    ) -> "RandomForestClassifierCustom":
        """
        Fits the Random Forest classifier to the training data.

        Args:
            X (np.ndarray): Training features.
            y (np.ndarray): Training labels.
            n_jobs (int): Number of jobs to run in parallel for both fit and predict.

        Returns:
            self: The instance of RandomForestClassifierCustom.
        """
        self.classes_ = np.unique(y)
        args = [
            (X, y, self.max_features, self.max_depth, self.random_state + i)
            for i in range(self.n_estimators)
        ]
        if n_jobs == 1:
            results = [self._fit_tree(arg) for arg in args]
        else:
            with Pool(n_jobs) as pool:
                results = pool.map(self._fit_tree, args)
        self.trees, self.feat_ids_by_tree = zip(*results)
        return self

    def _predict_tree_proba(
        self, args: Tuple[DecisionTreeClassifier, np.ndarray, np.ndarray]
    ) -> np.ndarray:
        """
        Predicts class probabilities for a given tree and feature set.

        Args:
            args (tuple): Contains the decision tree, feature indices, and the test dataset.

        Returns:
            np.ndarray: The predicted probabilities.
        """
        tree, feat_ids, X = args
        return tree.predict_proba(X[:, feat_ids])

    def predict_proba(self, X: np.ndarray, n_jobs: int = 1) -> np.ndarray:
        """
        Predicts class probabilities for each tree in the forest.

        Args:
            X (np.ndarray): Test features.
            n_jobs (int): Number of jobs to run in parallel.

        Returns:
            np.ndarray: The average predicted probabilities across all trees.
        """
        args = [
            (tree, self.feat_ids_by_tree[i], X) for i, tree in enumerate(self.trees)
        ]
        if n_jobs == 1:
            results = [self._predict_tree_proba(arg) for arg in args]
        else:
            with Pool(n_jobs) as pool:
                results = pool.map(self._predict_tree_proba, args)
        proba_sum = sum(results)
        averaged_proba = proba_sum / len(self.trees)
        return averaged_proba
