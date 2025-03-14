# Contributing to EirVS

Welcome, and thank you for considering contributing to EirVS! We’re excited to have your help. Whether you’re fixing bugs, improving documentation, or proposing new features, your contributions make a huge difference.

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

We welcome your contributions via pull requests (PRs). For more information on Pull requests, refer to: https://github.com/orgs/community/discussions/146509. Please follow the steps below:

### 1. Fork the Repository
1. Go to the EirVS GitHub repository.
2. Click the **Fork** button at the top-right corner to create a copy of the repository under your account.
![image](https://github.com/user-attachments/assets/7f9bc5a1-66b7-4bee-b9c2-697f7d15a37c)

### 2. Clone Your Fork

In your forked repository, clone the fork to your local machine, changing <your-username> in the link below to your Github account.
![image](https://github.com/user-attachments/assets/bea698ce-bdb8-4d1e-955b-95e0f92cae4a)

```bash
# Clone your fork to your local machine
git clone https://github.com/<your-username>/EirVS.git
cd EirVS
```

### 3. Create a New Branch
Create a new branch for your work. Use a descriptive name for your branch:
```bash
git checkout -b fix-issue-123
```

### 4. Make Changes

1. Make Your Changes
   Implement your proposed changes in the code or documentation. Make sure to follow the project's coding style and structure.  

2. Set Up Your Development Environment
   Ensure you have all the dependencies and environment configured.  

   ```bash
   # Create and activate the conda environment
   conda env create -f EirVS/environment.yml
   conda activate eirvs
   
   # Install the package in editable mode
   pip install -e EirVS
   ```

3. Run Tests Before Submitting
   Verify that all tests pass successfully to ensure your changes do not introduce regressions.  

   ```bash
   # Run tests
   python -m pytest test/
   ```

4. Update or Add Documentation  
   If your changes impact the functionality or introduce new features, make sure to update the relevant documentation.  

5. Run Tests Again After Changes  
   Re-run the tests to verify that your changes do not break any existing functionality. If your changes were intended to modify test results (e.g., by updating the `goldenData` for pytest), ensure these updates are included in your pull request.  

   ```bash
   # Run tests after making changes
   python -m pytest test/
   ```
---


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
1. Go to the original EirVS repository.
2. Click on the **Pull Requests** tab, then click **New Pull Request**.
   
![image](https://github.com/user-attachments/assets/88da4eb0-4013-47b7-96f9-4975728ad58e)

4. Choose your branch from your fork and compare it with the main branch of the original repository.

![image](https://github.com/user-attachments/assets/d65ac4b0-6b03-4459-b6d8-85baf9251489)

5. Provide the following details:
   - **Title**: A short description of the changes.
   - **Description**: A detailed explanation of your changes and the issue it addresses.
6. Choose at least two reviewers who you think are the most suitable for your changes.
7. Submit the pull request.

### 8. Review and Feedback
1. Your PR will require **at least 2 reviewers** to approve the changes before it can be merged.
2. Respond to feedback from reviewers promptly.
3. Make changes to your branch if needed and push updates to the same branch.

Once all requested changes are addressed, your PR will be merged by the maintainers.

---



Thank you for contributing to EirVS! Your efforts help make this project better for everyone. If you have any questions, feel free to reach out by opening an issue or contacting a maintainer.

