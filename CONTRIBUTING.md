# Contributing to MolSanitizer

Welcome, and thank you for considering contributing to MolSanitizer! We’re excited to have your help. Whether you’re fixing bugs, improving documentation, or proposing new features, your contributions make a huge difference.

---

## How to Open an Issue

If you find a bug, have a question, or want to suggest a new feature, please open an issue. Here’s how:

1. Navigate to the **Issues** tab in the GitHub repository.
2. Click on **New Issue**.
3. Choose an appropriate issue template (e.g., bug report, feature request).
4. Provide the following details:
   - **Title**: A clear and descriptive title.
   - **Description**: A detailed explanation of the issue, including the problem or feature idea.
   - **Steps to Reproduce**: (For bugs) Provide a minimal reproducible example or a step-by-step guide to reproduce the issue.
   - **Environment**: List relevant environment details, such as:
     - Operating system
     - Python version
     - RDKit version
     - Any other dependencies
   - **Expected Behavior**: (For bugs) What you expected to happen.
   - **Actual Behavior**: (For bugs) What actually happened.
   - **Screenshots/Logs**: Attach any relevant screenshots or error logs if applicable.
5. Submit the issue and wait for feedback from maintainers.

---

## Contributing via Forks and Pull Requests

We welcome your contributions via pull requests (PRs). Please follow the steps below:

### 1. Fork the Repository
1. Go to the MolSanitizer GitHub repository.
2. Click the **Fork** button at the top-right corner to create a copy of the repository under your account.

### 2. Clone Your Fork
```bash
# Clone your fork to your local machine
git clone https://github.com/<your-username>/MolSanitizer.git
cd MolSanitizer
```

### 3. Create a New Branch
Create a new branch for your work. Use a descriptive name for your branch:
```bash
git checkout -b fix-issue-123
```

### 4. Make Changes
1. Implement your changes in the code or documentation.
2. Run the tests to ensure everything works correctly.
   ```bash
   # Install dependencies
   pip install -r requirements.txt

   # Run tests
   pytest
   ```
3. Add or update documentation as needed.

### 5. Commit Your Changes
1. Stage your changes:
   ```bash
   git add .
   ```
2. Write a descriptive commit message:
   ```bash
   git commit -m "Fix issue #123: Correct handling of tautomeric states"
   ```

### 6. Push Your Changes
Push your branch to your forked repository:
```bash
git push origin fix-issue-123
```

### 7. Open a Pull Request
1. Go to the original MolSanitizer repository.
2. Click on the **Pull Requests** tab, then click **New Pull Request**.
3. Choose your branch from your fork and compare it with the main branch of the original repository.
4. Provide the following details:
   - **Title**: A short description of the changes.
   - **Description**: A detailed explanation of your changes and the issue it addresses.
5. Submit the pull request.

### 8. Review and Feedback
1. Your PR will require **at least 2 reviewers** to approve the changes before it can be merged.
2. Respond to feedback from reviewers promptly.
3. Make changes to your branch if needed and push updates to the same branch.

Once all requested changes are addressed, your PR will be merged by the maintainers.

---



Thank you for contributing to MolSanitizer! Your efforts help make this project better for everyone. If you have any questions, feel free to reach out by opening an issue or contacting a maintainer.

